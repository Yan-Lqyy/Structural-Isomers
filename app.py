#!/usr/bin/env python3
import time
import datetime
import logging
import json # Make sure json is imported
import sys
import os # Make sure os is imported for app.secret_key

from flask import Flask, request, render_template, url_for, Response, stream_with_context, jsonify

# Use the new core_logic package structure
from core_logic import (
    RDKIT_AVAILABLE, SMARTS_LIBRARY_DATA, get_smarts_categories,
    find_candidates_by_structure,
    filter_candidates_by_nmr,
    fetch_compound_details,
)
from typing import List, Optional, Generator, Callable, Any, Dict, Tuple

# --- Flask App Initialization ---
app = Flask(__name__)
app.secret_key = os.urandom(24) # Use a random secret key for session management

# --- CUSTOM JINJA FILTER ---
@app.template_filter('fromjson')
def fromjson_filter(value):
    """Deserializes a JSON string."""
    if value is None:
        return None
    try:
        return json.loads(value)
    except (json.JSONDecodeError, TypeError) as e:
        app.logger.error(f"Error decoding JSON in template filter: {e} for value: {value[:100]}...")
        return None # Or return an empty dict/list, or raise an error, depending on desired behavior

# --- Web Specific Constants (can be moved to a config file later) ---
MAX_CIDS_SUMMARY_LIST = 100  # Max CIDS to show in the "full list" on results page
INITIAL_DETAILS_TO_SHOW = 3  # Number of detailed results to show initially
LOAD_MORE_BATCH_SIZE = 3     # Number of details to load with "Load More" button
DETAIL_IMAGE_WIDTH = 250
DETAIL_IMAGE_HEIGHT = 250

# Get SMARTS categories using the function from core_logic
SMARTS_CATEGORIES = get_smarts_categories(SMARTS_LIBRARY_DATA)


@app.context_processor
def inject_global_vars():
    return dict(
        rdkit_available=RDKIT_AVAILABLE,
        now=datetime.datetime.now(datetime.timezone.utc), # Standard way for timezone-aware UTC
        app_version="1.3" # Example app version
    )

# --- Input Parsing Utilities ---
def parse_nmr_signals(signal_input: str) -> Optional[int]:
    if not signal_input: return None
    try:
        signals = int(signal_input)
        return signals if signals >= 0 else None
    except ValueError:
        return None

def parse_nmr_ratio(ratio_input: str) -> Optional[List[int]]:
    if not ratio_input: return None
    try:
        # Split, strip, convert to int, filter out non-positive, then sort
        parts = [int(p.strip()) for p in ratio_input.split(':') if p.strip()]
        valid_parts = [p for p in parts if p > 0]
        if len(valid_parts) != len(parts) or not valid_parts: # Check if all original parts were valid and list is not empty
            return None
        return sorted(valid_parts) # Store sorted for canonical comparison
    except ValueError:
        return None
    except Exception as e: # Catch any other parsing error
        app.logger.error(f"Unexpected error parsing NMR ratio '{ratio_input}': {e}")
        return None

# --- SSE Formatting ---
def format_sse(data: dict, event: Optional[str] = None) -> str:
    try:
        json_data = json.dumps(data)
    except TypeError as e:
        app.logger.error(f"Failed to serialize data for SSE: {data} - {e}")
        # Create an error payload to send to the client
        error_payload = {"type": "error", "message": f"Internal SSE serialization error: {e}"}
        json_data = json.dumps(error_payload)
    
    msg = f"data: {json_data}\n"
    if event:
        msg = f"event: {event}\n{msg}"
    return msg + "\n" # SSE messages must end with \n\n

# --- Routes ---
@app.route('/')
def index():
    app.logger.info("Rendering index page.")
    return render_template('index.html', categories=SMARTS_CATEGORIES)

@app.route('/search-stream')
def search_stream():
    formula = request.args.get('formula', '').strip()
    structural_constraint_keys = request.args.getlist('structural_constraints')
    nmr_signals_str = request.args.get('nmr_signals', '').strip()
    nmr_ratio_str = request.args.get('nmr_ratio', '').strip()

    @stream_with_context
    def generate_updates() -> Generator[str, None, None]:
        start_total_time = time.time()
        
        _messages_to_yield_immediately: List[str] = []
        _step_counter_for_client: int = 0

        def core_logic_progress_callback(message_from_core: str):
            nonlocal _step_counter_for_client
            _step_counter_for_client += 1
            payload = {"step": _step_counter_for_client, "type": "log", "message": message_from_core}
            try:
                _messages_to_yield_immediately.append(format_sse(payload))
            except Exception as e:
                app.logger.error(f"APP.PY (core_logic_callback): Error formatting/queuing progress: {e}")
                error_payload = {"step": _step_counter_for_client, "type": "error", "message": f"Internal server error processing progress: {e}"}
                _messages_to_yield_immediately.append(format_sse(error_payload))

        def queue_app_event(message: str, msg_type: str = "status", extra_data: Optional[dict] = None):
            nonlocal _step_counter_for_client
            _step_counter_for_client += 1
            payload = {"step": _step_counter_for_client, "type": msg_type, "message": message}
            if extra_data: payload.update(extra_data)
            try:
                _messages_to_yield_immediately.append(format_sse(payload))
            except Exception as e:
                app.logger.error(f"APP.PY (queue_app_event): Error formatting/queuing event: {e}")


        def yield_queued_messages():
            flushed_count = 0
            while _messages_to_yield_immediately:
                msg = _messages_to_yield_immediately.pop(0)
                yield msg
                flushed_count +=1
            if flushed_count > 0:
                sys.stdout.flush() # Ensure data is sent
                time.sleep(0.05) # Small pause after a burst of messages

        try:
            queue_app_event("Initiating search...", msg_type="start")
            yield from yield_queued_messages()

            # --- Input Validation ---
            queue_app_event("Validating input parameters...")
            validation_error_message = None
            client_warning_message = None # For warnings displayed on results page
            
            search_params_for_results = { # Store params for constructing results URL
                'formula': formula, 
                'structure_keys': ",".join(structural_constraint_keys) if structural_constraint_keys else "",
                'nmr_signals_str': nmr_signals_str, 
                'nmr_ratio_str': nmr_ratio_str,
                # Parsed versions will be added below if valid
            }

            if not formula: validation_error_message = "Molecular formula is required."
            elif not ('C' in formula.upper() and 'H' in formula.upper()): # Basic check
                validation_error_message = (validation_error_message + " " if validation_error_message else "") + \
                                       "Formula must contain at least Carbon (C) and Hydrogen (H)."
            
            nmr_signals = parse_nmr_signals(nmr_signals_str)
            if nmr_signals_str and nmr_signals is None:
                validation_error_message = (validation_error_message + " " if validation_error_message else "") + \
                                       "Invalid NMR signal count format (must be a non-negative integer)."
            search_params_for_results['nmr_signals'] = nmr_signals # Store parsed or None
            
            nmr_ratio = parse_nmr_ratio(nmr_ratio_str)
            if nmr_ratio_str and nmr_ratio is None:
                validation_error_message = (validation_error_message + " " if validation_error_message else "") + \
                                       "Invalid NMR ratio format (e.g., 1:2:3, integers > 0)."
            search_params_for_results['nmr_ratio'] = nmr_ratio # Store parsed or None

            if nmr_signals is not None and nmr_ratio is not None and len(nmr_ratio) != nmr_signals:
                 w_msg = "NMR signals count does not match the number of components in the ratio. Constraints will be applied independently if possible."
                 client_warning_message = w_msg
                 queue_app_event(w_msg, msg_type="warning")
            
            yield from yield_queued_messages()

            if validation_error_message:
                queue_app_event(validation_error_message, msg_type="error")
                # Prepare for completion event that will NOT redirect but allow user to see error
                completion_payload = {"type": "error_final", "message": validation_error_message}
                _messages_to_yield_immediately.append(format_sse(completion_payload, event="completion"))
                yield from yield_queued_messages()
                return

            # --- Prepare for Core Logic ---
            # Get actual SMARTS patterns for selected keys
            structural_constraints_smarts = {
                k: SMARTS_LIBRARY_DATA[k]['smarts'] 
                for k in structural_constraint_keys 
                if k in SMARTS_LIBRARY_DATA
            }

            # --- Execute Search: Step 1 (Formula + Structure) ---
            queue_app_event(f"Searching by formula '{formula}' and {len(structural_constraints_smarts)} structural constraint(s)...")
            yield from yield_queued_messages()

            candidate_cids, struct_search_error = find_candidates_by_structure(
                formula, structural_constraints_smarts, progress_callback=core_logic_progress_callback
            )
            yield from yield_queued_messages() # Yield messages from find_candidates_by_structure

            if struct_search_error:
                err_msg = f"Structural search failed: {struct_search_error}"
                queue_app_event(err_msg, msg_type="error")
                completion_payload = {"type": "error_final", "message": err_msg}
                _messages_to_yield_immediately.append(format_sse(completion_payload, event="completion"))
                yield from yield_queued_messages()
                return
            
            if candidate_cids is None: # Should be caught by struct_search_error, but defensive
                err_msg = "Unknown error during structure search, no CIDs returned."
                queue_app_event(err_msg, msg_type="error")
                completion_payload = {"type": "error_final", "message": err_msg}
                _messages_to_yield_immediately.append(format_sse(completion_payload, event="completion"))
                yield from yield_queued_messages()
                return

            if not candidate_cids:
                final_msg = "No compounds found matching the specified formula and structural criteria."
                queue_app_event(final_msg, msg_type="complete", extra_data={"count": 0})
                # Prepare for results page URL, even if no results
                results_url_params = search_params_for_results.copy()
                results_url_params['cids'] = "" # Empty CID list
                results_url_params['total_found'] = 0
                if client_warning_message: results_url_params['warning'] = client_warning_message
                
                final_url = url_for('results_display', **{k: v for k, v in results_url_params.items() if v is not None})
                completion_payload = {"type": "results_ready", "redirect_url": final_url, "message": final_msg, "count": 0}
                _messages_to_yield_immediately.append(format_sse(completion_payload, event="completion"))
                yield from yield_queued_messages()
                return

            # --- Execute Search: Step 2 (NMR Filtering, if applicable) ---
            final_cids_after_nmr = candidate_cids # Initialize
            nmr_filter_applied = False
            
            # Check if NMR filtering is actually requested and possible
            perform_nmr_filter = RDKIT_AVAILABLE and (nmr_signals is not None or nmr_ratio is not None)

            if perform_nmr_filter:
                nmr_filter_applied = True
                queue_app_event(f"Applying NMR filters to {len(candidate_cids)} candidates...")
                yield from yield_queued_messages()

                filtered_cids, nmr_filter_warning_msg = filter_candidates_by_nmr(
                    candidate_cids, nmr_signals, nmr_ratio, progress_callback=core_logic_progress_callback
                )
                yield from yield_queued_messages() # Yield messages from NMR filtering

                if nmr_filter_warning_msg:
                    # Prepend to existing client_warning_message
                    client_warning_message = (f"{nmr_filter_warning_msg}; {client_warning_message}" if client_warning_message else nmr_filter_warning_msg)
                    queue_app_event(f"NMR Filter Info: {nmr_filter_warning_msg}", msg_type="warning") # Also send as SSE event
                
                final_cids_after_nmr = filtered_cids
            else: # NMR filtering not performed
                if not RDKIT_AVAILABLE and (nmr_signals is not None or nmr_ratio is not None):
                     msg = "RDKit not available, NMR filters were skipped."
                     client_warning_message = (f"{msg}; {client_warning_message}" if client_warning_message else msg)
                     queue_app_event(msg, msg_type="warning")
                elif nmr_signals is None and nmr_ratio is None:
                      queue_app_event("NMR filtering skipped: No NMR constraints were specified.")
            
            yield from yield_queued_messages()


            # --- Final Check & Prepare Completion ---
            if not final_cids_after_nmr:
                 final_msg = "No compounds matched all criteria after NMR filtering." if nmr_filter_applied else "No compounds found matching structural criteria (NMR not applied or no effect)."
                 queue_app_event(final_msg, msg_type="complete", extra_data={"count": 0})
                 
                 results_url_params = search_params_for_results.copy()
                 results_url_params['cids'] = ""
                 results_url_params['total_found'] = 0
                 if client_warning_message: results_url_params['warning'] = client_warning_message
                 
                 final_url = url_for('results_display', **{k: v for k, v in results_url_params.items() if v is not None})
                 completion_payload = {"type": "results_ready", "redirect_url": final_url, "message": final_msg, "count": 0}
                 _messages_to_yield_immediately.append(format_sse(completion_payload, event="completion"))
                 yield from yield_queued_messages()
                 return

            # --- Success Path ---
            final_count = len(final_cids_after_nmr)
            success_msg = f"Search completed. Found {final_count} compound(s) matching all criteria."
            queue_app_event(success_msg, msg_type="complete", extra_data={"count": final_count})

            results_url_params = search_params_for_results.copy()
            
            # CIDs are passed to results page for initial display and "load more" logic
            # No need to truncate CIDs list passed to results_display anymore, as it handles display limits.
            results_url_params['cids'] = ",".join(map(str, final_cids_after_nmr))
            results_url_params['total_found'] = final_count # The actual number found
            
            if client_warning_message: results_url_params['warning'] = client_warning_message
            
            # Clean Nones before creating URL
            cleaned_results_params = {k: v for k, v in results_url_params.items() if v is not None}
            final_url_for_results = url_for('results_display', **cleaned_results_params)
            
            completion_payload = {
                "type": "results_ready", 
                "redirect_url": final_url_for_results, 
                "message": success_msg + " Click 'Show Results' to view.",
                "count": final_count
            }
            _messages_to_yield_immediately.append(format_sse(completion_payload, event="completion"))
            yield from yield_queued_messages()

            end_total_time = time.time()
            app.logger.info(f"SSE Stream processing completed in {end_total_time - start_total_time:.2f} seconds. Found {final_count} CIDs.")

        except Exception as e:
            app.logger.exception("Critical error during SSE stream generation (outer try-except):")
            # Attempt to yield a final error message to the client
            try:
                final_error_msg = f"An unexpected server error occurred: {str(e)}"
                # Use "error_final" type for client to know it's a terminal error from server side.
                completion_payload_error = {"type": "error_final", "message": final_error_msg}
                yield format_sse(completion_payload_error, event="completion")
                sys.stdout.flush()
            except Exception as e_final_yield:
                app.logger.error(f"Could not yield final error message to client after outer exception: {e_final_yield}")

    # Configure SSE response
    response = Response(generate_updates(), mimetype='text/event-stream')
    response.headers['Cache-Control'] = 'no-cache, no-transform' # Crucial for SSE
    response.headers['X-Accel-Buffering'] = 'no' # For Nginx, if used as a reverse proxy
    response.headers['Connection'] = 'keep-alive'
    return response


@app.route('/results')
def results_display():
    app.logger.info("Rendering results page.")
    # Retrieve search parameters from URL query
    formula = request.args.get('formula', '')
    structure_keys_str = request.args.get('structure_keys', '') # Comma-separated string
    nmr_signals_str_orig = request.args.get('nmr_signals_str', '') # Original input string
    nmr_ratio_str_orig = request.args.get('nmr_ratio_str', '')   # Original input string
    
    cids_str = request.args.get('cids', '') # Comma-separated string of all found CIDs
    total_found = request.args.get('total_found', type=int, default=0)
    
    warning_message_from_search = request.args.get('warning', None)
    error_message_from_search = request.args.get('error', None) # Should be rare if SSE handles errors

    # Parse parameters for display
    structure_keys_list = structure_keys_str.split(',') if structure_keys_str else []
    # Re-parse NMR constraints for consistent display (could also pass parsed versions if sure they are strings)
    nmr_signals_parsed_display = parse_nmr_signals(nmr_signals_str_orig)
    nmr_ratio_parsed_display = parse_nmr_ratio(nmr_ratio_str_orig)

    search_params_for_template = {
        'formula': formula,
        'structure_keys': structure_keys_list,
        'nmr_signals_str': nmr_signals_str_orig, # Show original input
        'nmr_ratio_str': nmr_ratio_str_orig,     # Show original input
        # For logic, we might prefer the parsed versions if they were passed directly
        # 'nmr_signals_parsed': nmr_signals_parsed_display,
        # 'nmr_ratio_parsed': nmr_ratio_parsed_display
    }

    all_found_cids: List[int] = []
    current_error = error_message_from_search # Prioritize errors from search process
    
    if cids_str:
        try:
            all_found_cids = [int(cid_item) for cid_item in cids_str.split(',') if cid_item.strip().isdigit()]
            # Validate if parsing matched the expected number based on split
            if len(all_found_cids) != len([c for c in cids_str.split(',') if c.strip()]):
                 parse_err = "Error parsing some CIDs from the list."
                 current_error = (current_error + "; " + parse_err) if current_error else parse_err
                 # Potentially clear all_found_cids or handle partially
        except ValueError:
            parse_err = "Invalid CID format in the received list."
            current_error = (current_error + "; " + parse_err) if current_error else parse_err
            all_found_cids = [] # Clear if major parsing error

    # If total_found from URL doesn't match parsed CIDs length (and no error yet), it might be an issue.
    if not current_error and total_found != len(all_found_cids):
        app.logger.warning(f"Mismatch: total_found ({total_found}) vs parsed CIDs ({len(all_found_cids)}). Using parsed CIDs length if available, else total_found.")
        # This could happen if CID string was truncated or malformed, but parsing succeeded partially.
        # total_found should be the ultimate source of truth for the "Found X" message.

    # Fetch initial batch of details
    detailed_results_initial: Dict[int, Dict[str, Any]] = {}
    cids_for_initial_detail = all_found_cids[:INITIAL_DETAILS_TO_SHOW]

    current_warning = warning_message_from_search

    if cids_for_initial_detail and not current_error:
        app.logger.info(f"Fetching initial details for display (up to {INITIAL_DETAILS_TO_SHOW} CIDs).")
        try:
            detailed_results_initial = fetch_compound_details(
                cids_for_initial_detail, DETAIL_IMAGE_WIDTH, DETAIL_IMAGE_HEIGHT
            )
            if len(detailed_results_initial) < len(cids_for_initial_detail):
                 detail_warn = f"Could only fetch details for {len(detailed_results_initial)} of the first {len(cids_for_initial_detail)} compounds."
                 current_warning = (current_warning + "; " + detail_warn) if current_warning else detail_warn
                 app.logger.warning(detail_warn)
        except Exception as e_detail:
             detail_err = f"Error fetching initial compound details: {e_detail}"
             current_error = (current_error + "; " + detail_err) if current_error else detail_err
             app.logger.exception("Error during initial detail fetching for results page:")
    
    # Ensure total_found is accurate for display message even if cids list is empty due to error
    display_total_found = total_found if total_found > 0 else len(all_found_cids)

    return render_template(
        'results.html',
        search_params=search_params_for_template,
        all_cids_json=json.dumps(all_found_cids), # Pass all CIDs as JSON for "load more" JS
        total_found=display_total_found,
        initial_detailed_results=detailed_results_initial,
        initial_display_count=len(detailed_results_initial),
        load_more_batch_size=LOAD_MORE_BATCH_SIZE,
        detail_image_width=DETAIL_IMAGE_WIDTH, # For JS to construct image URLs if needed
        detail_image_height=DETAIL_IMAGE_HEIGHT,
        max_cids_summary_list=MAX_CIDS_SUMMARY_LIST,
        error=current_error,
        warning=current_warning
    )

@app.route('/fetch-batch-details', methods=['POST'])
def fetch_batch_details_ajax():
    data = request.get_json()
    if not data or 'cids' not in data or not isinstance(data['cids'], list):
        return jsonify({"error": "Invalid request. 'cids' list is required."}), 400

    cids_to_fetch: List[int] = []
    try:
        cids_to_fetch = [int(cid) for cid in data['cids']]
    except ValueError:
        return jsonify({"error": "Invalid CIDs provided. Must be integers."}), 400
    
    if not cids_to_fetch:
        return jsonify({}) # Nothing to fetch

    app.logger.info(f"AJAX request to fetch details for {len(cids_to_fetch)} CIDs: {cids_to_fetch[:5]}...")
    try:
        details = fetch_compound_details(
            cids_to_fetch,
            request.args.get('img_width', DETAIL_IMAGE_WIDTH, type=int), # Allow overriding via query params if needed
            request.args.get('img_height', DETAIL_IMAGE_HEIGHT, type=int)
        )
        # The fetch_compound_details returns a dict {cid: details_dict}.
        # We might want to return it as a list of dicts in the order requested,
        # or ensure client side can handle dict. For simplicity, returning dict.

        # Prevent overloading server
        time.sleep(1) # Throttle requests to avoid overwhelming the server
        
        return jsonify(details)
    except Exception as e:
        app.logger.exception(f"Error in /fetch-batch-details for CIDs {cids_to_fetch[:5]}:")
        return jsonify({"error": f"Server error fetching details: {str(e)}"}), 500


# --- Main Execution ---
if __name__ == '__main__':
    print("-" * 60)
    print(" Starting PubChem Isomer Finder Web Application (v1.3)")
    print(" --> Using Gevent WSGI Server for robust SSE support <--")
    print("-" * 60)

    if not RDKIT_AVAILABLE:
        print("\n!!! WARNING: RDKit not found. NMR-related features will be disabled. !!!\n")
    else:
        print("\n+++ RDKit found and enabled. +++\n")

    # Configure logging
    log_level = logging.DEBUG if os.environ.get("FLASK_DEBUG") == "1" else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s %(levelname)s [%(name)s:%(lineno)d] %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        force=True # Override any existing handlers
    )
    
    # Set specific log levels for verbose libraries if needed
    logging.getLogger('werkzeug').setLevel(logging.INFO if log_level == logging.DEBUG else logging.WARNING)
    logging.getLogger('urllib3').setLevel(logging.INFO if log_level == logging.DEBUG else logging.WARNING)
    # Flask app's own logger will inherit from root or can be set explicitly
    app.logger.setLevel(log_level)
    # Core logic logger
    logging.getLogger('core_logic').setLevel(log_level)


    print(f"Flask App Log Level: {logging.getLevelName(app.logger.getEffectiveLevel())}")
    print(f"Core Logic Log Level: {logging.getLevelName(logging.getLogger('core_logic').getEffectiveLevel())}")
    print(f"Initial details to show: {INITIAL_DETAILS_TO_SHOW}, Load more batch: {LOAD_MORE_BATCH_SIZE}")


    flask_port = int(os.environ.get("PORT", 5001))
    try:
        from gevent.pywsgi import WSGIServer
        # app.debug = False # Gevent should run with debug False for production-like behavior
        print(f"Starting server with Gevent WSGIServer on http://0.0.0.0:{flask_port}")
        http_server = WSGIServer(('0.0.0.0', flask_port), app, log=logging.getLogger('gevent')) # Use a logger for gevent
        http_server.serve_forever()
    except ImportError:
        print("\n--- Gevent not found ---")
        print("Running with Flask's development server (SSE might be less reliable / buffered).")
        print("For optimal SSE, install gevent: pip install gevent")
        app.run(debug=(log_level == logging.DEBUG), host='0.0.0.0', port=flask_port, threaded=True) # threaded=True helps with multiple SSE clients
    except Exception as e:
        print(f"Failed to start server: {e}")
        sys.exit(1)