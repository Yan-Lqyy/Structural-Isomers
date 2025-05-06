#!/usr/bin/env python3
# app.py - Flask Web Application

import time
import datetime
import logging
from flask import Flask, request, render_template, url_for

# Import core logic functions and constants
from core_logic import (
    RDKIT_AVAILABLE, SMARTS_LIBRARY_DATA, get_smarts_categories,
    find_candidates_by_structure, filter_candidates_by_nmr,
    fetch_compound_details
)
from typing import List, Optional # Import necessary types

# --- Flask App Initialization ---
app = Flask(__name__)
# !!! CHANGE THIS IN PRODUCTION !!!
# Use environment variables or a config file for the secret key
app.secret_key = 'dev-secret-key-replace-me' # Keep simple for dev

# --- Web Specific Constants ---
MAX_CIDS_SUMMARY = 100
MAX_DETAILS_TO_SHOW = 10
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
    """Safely parse NMR signal count input from web form."""
    if not signal_input: return None
    try:
        signals = int(signal_input)
        return signals if signals >= 0 else None
    except ValueError:
        return None

def parse_nmr_ratio(ratio_input: str) -> Optional[List[int]]:
    """Safely parse NMR ratio input string into sorted list from web form."""
    if not ratio_input: return None
    try:
        parts = [int(p.strip()) for p in ratio_input.split(':') if p.strip()]
        valid_parts = [p for p in parts if p > 0]
        if len(valid_parts) != len(parts) or not valid_parts: return None
        return sorted(valid_parts)
    except ValueError: return None
    except Exception: return None

# --- Flask Routes ---

@app.route('/')
def index():
    """Renders the main search form."""
    app.logger.info("Rendering index page.")
    # Pass categories for building the form dynamically
    return render_template('index.html', categories=SMARTS_CATEGORIES)

@app.route('/search', methods=['POST'])
def search_results():
    """Handles form submission, performs search, renders results."""
    start_total_time = time.time()
    app.logger.info("Received search request.")

    # --- Get Form Data ---
    formula = request.form.get('formula', '').strip()
    structural_constraint_keys = request.form.getlist('structural_constraints')
    nmr_signals_str = request.form.get('nmr_signals', '').strip()
    nmr_ratio_str = request.form.get('nmr_ratio', '').strip()

    # --- Initialize Variables ---
    error_message = None
    warning_message = None
    final_cids = None
    detailed_results = {}
    search_params = {
        'formula': formula,
        'structure_keys': structural_constraint_keys,
        'nmr_signals': None,
        'nmr_ratio': None
    }

    # --- Input Validation ---
    if not formula:
        error_message = "Molecular formula is required."
    # Basic formula check (add more robust if needed)
    elif not ('C' in formula.upper() and 'H' in formula.upper()):
         error_message = "Formula must contain at least C and H."

    nmr_signals = parse_nmr_signals(nmr_signals_str)
    if nmr_signals_str and nmr_signals is None:
        error_message = (error_message + " " if error_message else "") + "Invalid NMR signal count entered."
    search_params['nmr_signals'] = nmr_signals

    nmr_ratio = parse_nmr_ratio(nmr_ratio_str)
    if nmr_ratio_str and nmr_ratio is None:
        error_message = (error_message + " " if error_message else "") + "Invalid NMR ratio format entered (use digits separated by colons, e.g., 3:1)."
    search_params['nmr_ratio'] = nmr_ratio

    # NMR Consistency Warning
    if nmr_signals is not None and nmr_ratio is not None and len(nmr_ratio) != nmr_signals:
         warning_message = "Warning: Number of NMR signals specified does not match the number of components in the ratio. Both constraints will be applied independently."

    # Map selected keywords to SMARTS strings for core logic
    structural_constraints_dict = {
        key: SMARTS_LIBRARY_DATA[key]['smarts']
        for key in structural_constraint_keys if key in SMARTS_LIBRARY_DATA
    }

    # --- Execute Search if Input is Valid ---
    if not error_message:
        app.logger.info(f"Executing search: Formula='{formula}', Structure Keys={structural_constraint_keys}, NMR Signals={nmr_signals}, NMR Ratio={nmr_ratio}")

        # 1. Find candidates by structure (call core logic)
        candidate_cids, struct_error = find_candidates_by_structure(formula, structural_constraints_dict)

        if struct_error:
            error_message = f"Structural search failed: {struct_error}"
        elif candidate_cids is None:
            error_message = "An unknown error occurred during the structural search phase."
        elif not candidate_cids:
            app.logger.info("No compounds found matching structural constraints.")
            final_cids = [] # Found zero, not an error
        else:
            # 2. Filter by NMR (call core logic)
            final_cids, nmr_warning = filter_candidates_by_nmr(candidate_cids, nmr_signals, nmr_ratio)
            if nmr_warning:
                warning_message = (warning_message + " " if warning_message else "") + nmr_warning

            # 3. Fetch details if results exist
            if final_cids:
                app.logger.info(f"Fetching details for {min(len(final_cids), MAX_DETAILS_TO_SHOW)} CIDs.")
                cids_for_detail = final_cids[:MAX_DETAILS_TO_SHOW]
                # Pass image dimensions needed by core logic now
                detailed_results = fetch_compound_details(
                    cids_for_detail,
                    DETAIL_IMAGE_WIDTH,
                    DETAIL_IMAGE_HEIGHT
                )

    # --- Render Results Template ---
    end_total_time = time.time()
    log_msg = f"Search request completed in {end_total_time - start_total_time:.2f} seconds."
    if error_message:
        app.logger.error(f"{log_msg} Error: {error_message}")
    else:
        app.logger.info(f"{log_msg} Found {len(final_cids) if final_cids is not None else 'N/A'} results.")

    return render_template(
        'results.html',
        search_params=search_params,
        final_cids=final_cids,
        detailed_results=detailed_results,
        max_cids_summary=MAX_CIDS_SUMMARY,
        error=error_message,
        warning=warning_message
    )

# --- Main Execution ---
if __name__ == '__main__':
    print("----------------------------------------------------")
    print(" Starting PubChem Structural Isomer Finder Web App")
    print("----------------------------------------------------")
    if not RDKIT_AVAILABLE:
        print("\n!!! WARNING: RDKit not found. NMR features will be disabled. !!!\n")
    # Enable Flask logging for development
    app.logger.setLevel(logging.INFO)
    # Run in debug mode for development (auto-reload, debugger)
    # Use host='0.0.0.0' carefully to make accessible on local network
    app.run(debug=True, host='0.0.0.0', port=5001) # Use a different port if 5000 is busy
    # For production deployment, use a proper WSGI server like Gunicorn:
    # gunicorn -w 4 -b 0.0.0.0:5001 app:app