import logging
from typing import List, Dict, Optional, Tuple, Callable
import urllib.parse

from .pubchem_api import run_pubchem_search, BASE_URL

logger = logging.getLogger(__name__)

def find_candidates_by_structure(
    formula: str,
    constraints_dict: Dict[str, str],
    progress_callback: Optional[Callable[[str], None]] = None
) -> Tuple[Optional[List[int]], Optional[str]]:
    """
    Finds candidate CIDs by formula and structural constraints using PubChem PUG REST.
    Returns a tuple: (list_of_cids, error_message_or_none).
    """
    def report_progress(message: str):
        logger.info(f"[StructSearch] {message}")
        if progress_callback:
            try:
                progress_callback(message)
            except Exception as e:
                logger.warning(f"Progress callback failed for message '{message[:50]}...': {e}")

    report_progress(f"Initiating structure search for formula: '{formula}'")

    # Step 1: Search by formula to get an initial list (or cachekey)
    formula_url = f"{BASE_URL}/compound/fastformula/{formula}/cids/JSON"
    # Request cachekey to allow subsequent refinements
    formula_params = {'list_return': 'cachekey'} 
    formula_data = run_pubchem_search(formula_url, params=formula_params)

    if '_error_status' in formula_data:
        msg = formula_data.get('_error_message', 'Unknown error during formula search')
        report_progress(f"Formula search failed: {msg}")
        return None, msg
    
    if not formula_data or 'IdentifierList' not in formula_data or 'CacheKey' not in formula_data['IdentifierList']:
        # This case could also mean 0 results if PubChem doesn't return 'Size' for empty cachekeys
        size_info = formula_data.get('IdentifierList', {}).get('Size', 'N/A')
        if size_info == 0 :
            report_progress(f"No compounds found matching formula '{formula}'.")
            return [], None # No error, just no results
        msg = "Invalid response or missing CacheKey from formula search."
        report_progress(msg)
        return None, msg

    current_cachekey = formula_data['IdentifierList']['CacheKey']
    initial_size = formula_data['IdentifierList'].get('Size', 0) # Size might not be present if 0
    
    if initial_size == 0:
        report_progress(f"No compounds found matching formula '{formula}'.")
        return [], None

    report_progress(f"Formula '{formula}' matched {initial_size} compound(s). (CacheKey: {current_cachekey[:15]}...)")

    # Step 2: Apply structural constraints if any
    if not constraints_dict:
        report_progress("No structural constraints. Retrieving all CIDs from formula search cache...")
        # Retrieve all CIDs from the cachekey
        final_cids_url = f"{BASE_URL}/compound/listkey/{current_cachekey}/cids/JSON"
        final_cids_data = run_pubchem_search(final_cids_url)

        if '_error_status' in final_cids_data:
            msg = final_cids_data.get('_error_message', 'Error retrieving full CIDs list from cache')
            report_progress(f"Failed to retrieve CIDs: {msg}")
            return None, msg
        if final_cids_data and 'IdentifierList' in final_cids_data and 'CID' in final_cids_data['IdentifierList']:
            cids = final_cids_data['IdentifierList']['CID']
            report_progress(f"Retrieved {len(cids)} CIDs matching formula only.")
            return cids, None
        else:
            msg = "Error parsing CIDs from formula cachekey results."
            report_progress(msg)
            return None, msg

    # Apply constraints one by one
    report_progress(f"Applying {len(constraints_dict)} structural constraint(s) sequentially...")
    num_constraints = len(constraints_dict)
    constraint_items = list(constraints_dict.items()) 

    for i, (constraint_key, smarts_pattern) in enumerate(constraint_items):
        is_last_constraint = (i == num_constraints - 1)
        
        # --- URL ENCODE THE SMARTS PATTERN ---
        encoded_smarts_pattern = urllib.parse.quote(smarts_pattern)
        # --------------------------------------

        report_progress(f"  Applying constraint {i+1}/{num_constraints}: '{constraint_key}' (Original SMARTS: {smarts_pattern[:30]}...)")
        
        # Use the encoded SMARTS pattern in the URL
        refine_url = f"{BASE_URL}/compound/fastsubstructure/smarts/{encoded_smarts_pattern}/cids/JSON"
        refine_params = {'cachekey': current_cachekey}
        
        expected_response_key = 'CID' if is_last_constraint else 'CacheKey'
        if not is_last_constraint:
            refine_params['list_return'] = 'cachekey'

        refine_data = run_pubchem_search(refine_url, params=refine_params)
        
        if '_error_status' in refine_data:
            msg = refine_data.get('_error_message', f"Error during SMARTS refinement for '{constraint_key}'")
            report_progress(f"Constraint application failed for '{constraint_key}': {msg}")
            return None, msg

        if not refine_data or 'IdentifierList' not in refine_data:
            msg = f"Unexpected empty or malformed response during refinement for '{constraint_key}'."
            report_progress(msg)
            return None, msg
        
        identifier_list = refine_data['IdentifierList']
        result_value = identifier_list.get(expected_response_key)
        current_size = identifier_list.get('Size', 0)

        if result_value is None and current_size == 0: # Constraint yielded 0 results
            report_progress(f"  Constraint '{constraint_key}' resulted in 0 matches. No compounds meet all criteria.")
            return [], None # No error, just no results
        
        if result_value is None: # Should not happen if size > 0
            msg = f"Missing '{expected_response_key}' in response for '{constraint_key}', though size is {current_size}."
            report_progress(msg)
            return None, msg

        if not is_last_constraint:
            current_cachekey = result_value
            report_progress(f"  Constraint '{constraint_key}' applied. {current_size} compound(s) remaining. (New CacheKey: {current_cachekey[:15]}...)")
        else:
            # This is the final list of CIDs after all structural constraints
            final_cids = result_value 
            report_progress(f"  Final structural constraint '{constraint_key}' applied. Found {len(final_cids)} candidate CID(s).")
            return final_cids, None
            
    # This part should ideally not be reached if logic is correct
    report_progress("Constraint application loop completed unexpectedly.")
    return None, "Internal error in structural constraint application logic."