#!/usr/bin/env python3
# app.py - Flask Web Application

import time
import datetime
import logging
import json
import sys # *** Import sys for explicit flushing ***
# *** Import necessary Flask components for SSE and streaming ***
from flask import Flask, request, render_template, url_for, Response, stream_with_context, jsonify

# *** Import core logic functions and constants, including updated ones for callback ***
# Ensure your core_logic.py has imported 'Callable' from 'typing'
# and moved NMR helper function definitions before filter_candidates_by_nmr
from core_logic import (
    RDKIT_AVAILABLE, SMARTS_LIBRARY_DATA, get_smarts_categories,
    find_candidates_by_structure, # This should now accept progress_callback
    filter_candidates_by_nmr,     # This should now accept progress_callback
    fetch_compound_details        # This does NOT need callback as it runs after the main search
)
# *** Import typing for the callback if needed in app.py directly (though usually only needed in core_logic params) ***
from typing import List, Optional, Generator, Callable, Any, Dict, Tuple

# --- Flask App Initialization ---
app = Flask(__name__)
# !!! CHANGE THIS IN PRODUCTION !!!
# Use environment variables or a config file for the secret key
app.secret_key = 'dev-secret-key-replace-me' # Keep simple for dev

# --- Web Specific Constants ---
MAX_CIDS_SUMMARY = 100
MAX_DETAILS_TO_SHOW = 5 # Keep low to avoid image rate limits
DETAIL_IMAGE_WIDTH = 250
DETAIL_IMAGE_HEIGHT = 250

# --- Get Categorized SMARTS for Form ---
# Compute this once at startup
SMARTS_CATEGORIES = get_smarts_categories(SMARTS_LIBRARY_DATA)

# --- Context Processor ---
@app.context_processor
def inject_global_vars():
    """Make selected variables available to all templates."""
    return dict(
        rdkit_available=RDKIT_AVAILABLE,
        now=datetime.datetime.utcnow()
    )

# --- Input Parsing Helpers ---
def parse_nmr_signals(signal_input: str) -> Optional[int]:
    if not signal_input: return None
    try:
        signals = int(signal_input)
        return signals if signals >= 0 else None
    except ValueError: return None

def parse_nmr_ratio(ratio_input: str) -> Optional[List[int]]:
    if not ratio_input: return None
    try:
        parts = [int(p.strip()) for p in ratio_input.split(':') if p.strip()]
        valid_parts = [p for p in parts if p > 0]
        if len(valid_parts) != len(parts) or not valid_parts: return None
        return sorted(valid_parts)
    except ValueError: return None
    except Exception: return None

# --- Helper to format SSE messages ---
def format_sse(data: dict, event: Optional[str] = None) -> str:
    """Formats data as a Server-Sent Event string."""
    try:
        json_data = json.dumps(data)
    except TypeError as e:
        app.logger.error(f"Failed to serialize data for SSE: {data} - {e}")
        json_data = json.dumps({"type": "error", "message": f"Serialization error: {e}"})
    msg = f"data: {json_data}\n"
    if event: msg = f"event: {event}\n{msg}"
    return msg + "\n"

# --- Flask Routes ---

@app.route('/')
def index():
    """Renders the main search form."""
    app.logger.info("Rendering index page.")
    return render_template('index.html', categories=SMARTS_CATEGORIES)

# --- SSE Search Stream Route ---
@app.route('/search-stream')
def search_stream():
    """Handles search via SSE, yielding progress updates."""
    formula = request.args.get('formula', '').strip()
    structural_constraint_keys = request.args.getlist('structural_constraints')
    nmr_signals_str = request.args.get('nmr_signals', '').strip()
    nmr_ratio_str = request.args.get('nmr_ratio', '').strip()

    @stream_with_context
    def generate_updates() -> Generator[str, None, None]:
        start_total_time = time.time()
        step_counter_for_client = 0

        # --- Define the progress callback function ---
        def core_logic_progress_callback(message: str):
            nonlocal step_counter_for_client
            step_counter_for_client += 1
            payload = {"step": step_counter_for_client, "type": "log", "message": message}
            app.logger.debug(f"Callback->SSE Yield: {message}") # Log the message being yielded
            yield format_sse(payload)
            sys.stdout.flush() # Force flush
            # Minimal sleep, might prevent overwhelming the client/network buffer
            time.sleep(0.02)

        # --- Define yield functions for stream-level events ---
        def yield_stream_event(message: str, msg_type: str = "status", extra_data: Optional[dict] = None):
             nonlocal step_counter_for_client
             step_counter_for_client += 1
             payload = {"step": step_counter_for_client, "type": msg_type, "message": message}
             if extra_data: payload.update(extra_data)
             app.logger.info(f"StreamEvent->SSE Yield ({msg_type}): {message}") # Log the message being yielded
             yield format_sse(payload)
             sys.stdout.flush() # Force flush
             time.sleep(0.05)

        def yield_error(message: str):
             app.logger.error(f"SSE Error Yield: {message}")
             yield format_sse({"type": "error", "message": message}, event="error")
             sys.stdout.flush() # Force flush
             time.sleep(0.05)

        # --- Start Actual Logic ---
        try:
            yield_stream_event("Initiating search...", msg_type="start")

            # --- Input Validation ---
            yield_stream_event("Validating input parameters...")
            error_message = None
            warning_message = None
            search_params = {
                'formula': formula, 'structure_keys': structural_constraint_keys,
                'nmr_signals': None, 'nmr_ratio': None,
                'nmr_signals_str': nmr_signals_str, 'nmr_ratio_str': nmr_ratio_str
            }
            # ... (Validation logic - same as before) ...
            if not formula: error_message = "Molecular formula is required."
            elif not ('C' in formula.upper() and 'H' in formula.upper()): error_message = "Formula must contain at least C and H."
            nmr_signals = parse_nmr_signals(nmr_signals_str)
            if nmr_signals_str and nmr_signals is None: error_message = (error_message + " " if error_message else "") + "Invalid NMR signal count entered."
            search_params['nmr_signals'] = nmr_signals
            nmr_ratio = parse_nmr_ratio(nmr_ratio_str)
            if nmr_ratio_str and nmr_ratio is None: error_message = (error_message + " " if error_message else "") + "Invalid NMR ratio format entered."
            search_params['nmr_ratio'] = nmr_ratio
            if nmr_signals is not None and nmr_ratio is not None and len(nmr_ratio) != nmr_signals:
                 w_msg = "Warning: NMR signals count doesn't match ratio components. Applied independently."
                 warning_message = w_msg; yield_stream_event(warning_message, msg_type="warning")
            if error_message: yield_error(error_message); return

            structural_constraints_dict = {k: SMARTS_LIBRARY_DATA[k]['smarts'] for k in structural_constraint_keys if k in SMARTS_LIBRARY_DATA}

            # --- Execute Search: Step 1 (Formula + Structure) ---
            yield_stream_event("Searching by formula and structure...") # Frame the step
            candidate_cids, struct_error = find_candidates_by_structure(
                formula, structural_constraints_dict, progress_callback=core_logic_progress_callback
            )
            if struct_error: yield_error(f"Structural search failed: {struct_error}"); return
            if candidate_cids is None: yield_error("Unknown error during structure search."); return

            if not candidate_cids:
                final_msg = "Search complete: No compounds found matching structural criteria."
                yield_stream_event(final_msg, msg_type="complete", extra_data={"count": 0})
                # ... prepare URL and yield completion redirect (same as before) ...
                results_params = search_params.copy(); results_params['cids'] = ""; results_params['total_found'] = 0; results_params['warning'] = warning_message
                final_url = url_for('results_display', **{k: v for k, v in results_params.items() if v is not None})
                yield format_sse({"type": "redirect", "url": final_url}, event="completion")
                return

            # --- Execute Search: Step 2 (NMR Filtering) ---
            final_cids = candidate_cids
            nmr_applied = False
            if RDKIT_AVAILABLE and (search_params['nmr_signals'] is not None or search_params['nmr_ratio'] is not None):
                nmr_applied = True
                yield_stream_event("Applying NMR filters...") # Frame the step
                filtered_cids, nmr_warning = filter_candidates_by_nmr(
                    candidate_cids, search_params['nmr_signals'], search_params['nmr_ratio'], progress_callback=core_logic_progress_callback
                )
                if nmr_warning:
                    full_nmr_warning = (warning_message + "; " if warning_message else "") + nmr_warning
                    warning_message = full_nmr_warning
                    yield_stream_event(f"NMR Warning: {nmr_warning}", msg_type="warning")
                final_cids = filtered_cids
            else:
                # ... logic for skipping NMR (same as before) ...
                if not RDKIT_AVAILABLE and (search_params['nmr_signals'] is not None or search_params['nmr_ratio'] is not None):
                     msg = "RDKit not available, NMR filters skipped."; warning_message = (warning_message + "; " if warning_message else "") + msg
                     yield_stream_event(msg, msg_type="warning")
                elif search_params['nmr_signals'] is None and search_params['nmr_ratio'] is None:
                      yield_stream_event("NMR filtering skipped: No constraints specified.")

            # --- Final Check & Prepare Completion ---
            if not final_cids:
                 final_msg = "Search complete: No compounds matched all criteria after filtering." if nmr_applied else "Search complete: No compounds found matching structural criteria."
                 yield_stream_event(final_msg, msg_type="complete", extra_data={"count": 0})
                 # ... prepare URL and yield completion redirect (same as before) ...
                 results_params = search_params.copy(); results_params['cids'] = ""; results_params['total_found'] = 0; results_params['warning'] = warning_message
                 final_url = url_for('results_display', **{k: v for k, v in results_params.items() if v is not None})
                 yield format_sse({"type": "redirect", "url": final_url}, event="completion")
                 return

            # Success path
            final_count = len(final_cids)
            yield_stream_event(f"Search completed successfully. Found {final_count} compounds.", msg_type="complete", extra_data={"count": final_count})

            # --- Generate URL for results page (same logic as before) ---
            results_params = search_params.copy()
            cids_limit_url = 500
            cids_to_pass = final_cids[:cids_limit_url]
            if len(final_cids) > len(cids_to_pass):
                 url_trunc_msg = f"Result list truncated in URL (passing {len(cids_to_pass)} of {len(final_cids)} CIDs)."
                 warning_message = (warning_message + "; " if warning_message else "") + url_trunc_msg
            results_params['cids'] = ",".join(map(str, cids_to_pass))
            results_params['total_found'] = len(final_cids)
            results_params['warning'] = warning_message
            results_params_cleaned = {k: v for k, v in results_params.items() if v is not None}
            final_url = url_for('results_display', **results_params_cleaned)
            yield format_sse({"type": "redirect", "url": final_url}, event="completion")

            end_total_time = time.time()
            app.logger.info(f"SSE Stream processing completed in {end_total_time - start_total_time:.2f} seconds.")

        except Exception as e:
            error_msg = f"An unexpected server error occurred during search: {e}"
            yield_error(error_msg)
            app.logger.exception("Error during SSE stream generation:")

    # Return the generator function wrapped in a Response object
    response = Response(generate_updates(), mimetype='text/event-stream')
    response.headers['Cache-Control'] = 'no-cache, no-transform'
    response.headers['X-Accel-Buffering'] = 'no'
    response.headers['Connection'] = 'keep-alive'
    return response


# --- GET Route to Display Results ---
# (No changes needed in this route, its logic is fine)
@app.route('/results')
def results_display():
    """Displays the final search results page."""
    app.logger.info("Rendering results page.")
    formula = request.args.get('formula', '')
    structure_keys_str = request.args.get('structure_keys', '')
    nmr_signals_str = request.args.get('nmr_signals_str', '')
    nmr_ratio_str = request.args.get('nmr_ratio_str', '')
    cids_str = request.args.get('cids', '')
    total_found = request.args.get('total_found', type=int, default=0)
    warning_message = request.args.get('warning', None)
    error_message = request.args.get('error', None)

    structure_keys = structure_keys_str.split(',') if structure_keys_str else []
    nmr_signals_parsed = parse_nmr_signals(nmr_signals_str)
    nmr_ratio_parsed = parse_nmr_ratio(nmr_ratio_str)
    search_params_display = {
        'formula': formula, 'structure_keys': structure_keys,
        'nmr_signals': nmr_signals_parsed, 'nmr_ratio': nmr_ratio_parsed,
        'nmr_signals_str': nmr_signals_str, 'nmr_ratio_str': nmr_ratio_str
    }
    final_cids = []
    if cids_str:
        try:
            final_cids = [int(cid) for cid in cids_str.split(',') if cid.isdigit()]
            if len(final_cids) != len([c for c in cids_str.split(',') if c]):
                 error_message = (error_message + " " if error_message else "") + "Error parsing CID list."
                 final_cids = []
        except ValueError:
            error_message = (error_message + " " if error_message else "") + "Invalid CID format received."
            final_cids = []
    detailed_results = {}
    if final_cids and not error_message:
        app.logger.info(f"Fetching details for display (top {min(len(final_cids), MAX_DETAILS_TO_SHOW)} CIDs).")
        cids_for_detail = final_cids[:MAX_DETAILS_TO_SHOW]
        try:
            detailed_results = fetch_compound_details(
                cids_for_detail, DETAIL_IMAGE_WIDTH, DETAIL_IMAGE_HEIGHT
            )
            if len(detailed_results) < len(cids_for_detail) and cids_for_detail:
                 det_warning = f"Could not fetch details for {len(cids_for_detail) - len(detailed_results)} top compounds."
                 warning_message = (warning_message + "; " if warning_message else "") + det_warning
                 app.logger.warning(det_warning)
        except Exception as e:
             error_message = (error_message + " " if error_message else "") + f"Error fetching details: {e}"
             app.logger.exception("Error during detail fetching:")
    if total_found != len(final_cids) and final_cids:
         if warning_message is None or "Result list truncated in URL" not in warning_message:
             trunc_warning = f"{total_found} results found, list in URL truncated to {len(final_cids)} CIDs."
             warning_message = (warning_message + "; " if warning_message else "") + trunc_warning
    if not final_cids and total_found > 0:
        warning_message = (warning_message + " " if warning_message else "") + f"{total_found} results found, list too long to pass/display."

    return render_template(
        'results.html', search_params=search_params_display, final_cids=final_cids,
        total_found=total_found, detailed_results=detailed_results,
        max_cids_summary=MAX_CIDS_SUMMARY, error=error_message, warning=warning_message
    )

# --- Main Execution (Using Gevent for better SSE handling) ---
if __name__ == '__main__':
    print("-" * 50)
    print(" Starting PubChem Structural Isomer Finder Web App")
    print(" --> Using Gevent WSGI Server for SSE <--")
    print("-" * 50)

    if not RDKIT_AVAILABLE:
        print("\n!!! WARNING: RDKit not found. NMR features will be disabled. !!!\n")

    # Configure logging (DEBUG for development, INFO for production)
    # Ensure logs are visible
    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(levelname)s [%(name)s:%(lineno)d] %(message)s', force=True)
    app.logger.setLevel(logging.DEBUG)
    logging.getLogger('werkzeug').setLevel(logging.INFO) # Werkzeug logger won't be used, but set anyway
    logging.getLogger('core_logic').setLevel(logging.DEBUG)
    # Gevent logger is usually quiet, but can be configured if needed
    # logging.getLogger('gevent').setLevel(logging.DEBUG)

    print(f"Flask App Log Level: {logging.getLevelName(app.logger.level)}")
    print(f"Core Logic Log Level: {logging.getLevelName(logging.getLogger('core_logic').level)}")
    print(f"Max details to show (images): {MAX_DETAILS_TO_SHOW}")

    # *** Use Gevent WSGI Server ***
    try:
        from gevent.pywsgi import WSGIServer
        # Note: Flask's 'debug=True' auto-reloader doesn't work well with Gevent.
        # You'll need to manually restart the server after code changes.
        # Set Flask debug=False when using WSGIServer, but keep logger level DEBUG.
        app.debug = False
        print("Starting server with Gevent WSGIServer on http://0.0.0.0:5001")
        http_server = WSGIServer(('0.0.0.0', 5001), app, log=app.logger) # Pass Flask logger to Gevent
        http_server.serve_forever()
    except ImportError:
        print("\n--- Gevent not found ---")
        print("Running with Flask's development server (SSE might be buffered).")
        print("For better SSE support, install gevent: pip install gevent")
        # Fallback to Werkzeug if Gevent is not installed
        app.run(debug=True, host='0.0.0.0', port=5001, threaded=True)