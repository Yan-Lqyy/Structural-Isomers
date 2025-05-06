import logging
import math
from collections import Counter
from typing import List, Dict, Optional, Tuple, Callable, Any # Added Any

from . import RDKIT_AVAILABLE # Use __init__.py for RDKIT_AVAILABLE
if RDKIT_AVAILABLE:
    from rdkit import Chem
    # from rdkit.Chem import AllChem # Already in __init__

from .compound_utils import fetch_smiles_batch # Import from sibling module

logger = logging.getLogger(__name__)

# --- NMR Prediction Helpers (RDKit dependent) ---
# --- NMR Prediction Helpers (RDKit dependent) ---
def get_hydrogen_environments(smiles_string: str) -> Tuple[Optional[int], Optional[Dict[int, int]]]:
    """
    Predicts the number of unique hydrogen environments (signals) and their integrations.
    Returns: (number_of_signals, {rank_id: count_of_hydrogens_for_rank}) or (None, None) on error.
    """
    if not RDKIT_AVAILABLE:
        logger.debug("RDKit not available, cannot get hydrogen environments.")
        return None, None
    if not smiles_string:
        logger.warning("get_hydrogen_environments called with empty SMILES string.")
        return None, None

    try:
        mol = Chem.MolFromSmiles(smiles_string)
        if mol is None:
            logger.debug(f"RDKit failed to parse SMILES: '{smiles_string[:50]}...'")
            return None, None

        mol_h = Chem.AddHs(mol)
        if mol_h is None: 
            logger.warning(f"RDKit failed to add Hydrogens to molecule from SMILES: '{smiles_string[:50]}...'")
            return None, None
        
        has_hydrogens = any(atom.GetAtomicNum() == 1 for atom in mol_h.GetAtoms())
        if not has_hydrogens:
            logger.debug(f"No hydrogens found in molecule from SMILES: '{smiles_string[:50]}...'")
            return 0, {} 

        # Get canonical ranks for all atoms (including hydrogens)
        # *** CORRECTED PARAMETER NAMES HERE ***
        ranks = list(Chem.CanonicalRankAtoms(mol_h, 
                                             breakTies=False, 
                                             includeChirality=True, # Changed from doChirality
                                             includeIsotopes=False  # Changed from doIsotopes
                                             # other params will use defaults: includeAtomMaps=True, includeChiralPresence=False
                                             ))
        
        h_counts_by_rank = Counter()
        for atom in mol_h.GetAtoms():
            if atom.GetAtomicNum() == 1:  
                try:
                    rank_id = ranks[atom.GetIdx()]
                    h_counts_by_rank[rank_id] += 1
                except IndexError:
                    logger.error(f"IndexError accessing ranks for atom {atom.GetIdx()} in SMILES: {smiles_string}")
                    continue 

        num_signals = len(h_counts_by_rank)
        integration_counts_dict = dict(h_counts_by_rank)
        
        if has_hydrogens and num_signals == 0 and h_counts_by_rank: # Typo: should be `not h_counts_by_rank` or just check `num_signals == 0`
             logger.warning(f"Molecule '{smiles_string[:50]}...' has hydrogens, but RDKit reported 0 unique environments. This might indicate an issue or a highly symmetrical molecule with no distinct H signals by this method.")
        
        logger.debug(f"NMR Prediction for '{smiles_string[:30]}...': {num_signals} signal(s), Integrations: {integration_counts_dict}")
        return num_signals, integration_counts_dict

    except Exception as e:
        # Log the full traceback for these C++ signature mismatch errors to help debug
        logger.warning(f"RDKit error processing SMILES '{smiles_string[:50]}...' for H environments: {e}", exc_info=True) 
        return None, None

def _gcd(a: int, b: int) -> int:
    """Helper for Greatest Common Divisor."""
    while(b):
        a, b = b, a % b
    return abs(a) # Ensure positive GCD

def _gcd_list(numbers: List[int]) -> int:
    """Calculates GCD of a list of numbers."""
    if not numbers:
        return 1 # Or raise error, but 1 is neutral for division
    
    # Filter out zeros or negative numbers if they shouldn't be part of ratio calc
    positive_numbers = [n for n in numbers if n > 0]
    if not positive_numbers:
        return 1 # All numbers were <=0

    result = positive_numbers[0]
    for i in range(1, len(positive_numbers)):
        result = _gcd(result, positive_numbers[i])
    return result if result > 0 else 1 # Ensure GCD is at least 1

def calculate_predicted_ratio(integration_counts_dict: Optional[Dict[Any, int]]) -> Optional[List[int]]:
    """
    Calculates the simplified proton ratio from integration counts.
    Returns a sorted list of integers representing the ratio, or None on error.
    """
    if integration_counts_dict is None:
        return None
    if not integration_counts_dict: # Empty dict means 0 signals, so empty ratio
        return []

    counts = list(integration_counts_dict.values())
    if not counts: # Should be caught by empty dict, but good check
        return []

    # Ensure all counts are positive for ratio calculation
    positive_counts = [c for c in counts if c > 0]
    if not positive_counts:
        logger.debug(f"Ratio calculation: No positive integration counts found in {counts}. Returning empty ratio.")
        return []

    common_divisor = _gcd_list(positive_counts)
    if common_divisor <= 0: # Should be handled by _gcd_list returning 1
        logger.error(f"Invalid GCD ({common_divisor}) calculated for counts: {positive_counts}. Cannot determine ratio.")
        return None 

    # Divide each count by the GCD and sort for a canonical representation
    ratio_list = sorted([count // common_divisor for count in positive_counts])
    
    logger.debug(f"Calculated ratio {ratio_list} from integration counts {counts} (GCD: {common_divisor})")
    return ratio_list


def format_nmr_prediction(num_signals: Optional[int], integration_counts_dict: Optional[Dict[Any, int]]) -> str:
    """Formats NMR prediction data into a human-readable string."""
    if num_signals is None:
        return "N/A (RDKit Prediction Error)"
    if num_signals == 0:
        return "0 Signals (No Chemically Distinct Hydrogens Expected or No Hydrogens)"

    integration_str = "N/A"
    if integration_counts_dict:
        # Display raw integration counts, sorted for consistency
        counts_display = sorted(list(integration_counts_dict.values()), reverse=True)
        integration_str = ", ".join([f"{count}H" for count in counts_display if count > 0])
    else: # num_signals > 0 but no integration_counts_dict (should not happen if logic is correct)
        integration_str = "Counts unavailable"


    ratio_list = calculate_predicted_ratio(integration_counts_dict)
    ratio_str = "N/A"
    if ratio_list is None: # Error in calculation
        ratio_str = "Error"
    elif not ratio_list and num_signals > 0 : # Non-zero signals but empty ratio (e.g. all counts were zero)
        ratio_str = "Undefined (check integrations)"
    elif ratio_list:
        ratio_str = ":".join(map(str, ratio_list))
    else: # num_signals == 0 implies empty ratio_list, this is fine
        ratio_str = "N/A (0 signals)"


    return f"{num_signals} Signal(s); Ratio: {ratio_str}; Integrations: ({integration_str})"


# --- NMR Filtering Logic ---
def filter_candidates_by_nmr(
    candidate_cids: List[int],
    required_signals: Optional[int],
    required_ratio: Optional[List[int]], # Assumed to be sorted
    progress_callback: Optional[Callable[[str], None]] = None
) -> Tuple[List[int], Optional[str]]: # Returns (filtered_cids, warning_message_or_none)
    """
    Filters a list of CIDs based on predicted NMR signals and ratio.
    """
    def report_progress(message: str):
        logger.info(f"[NMRFilter] {message}")
        if progress_callback:
            try:
                progress_callback(message)
            except Exception as e:
                logger.warning(f"Progress callback failed for message '{message[:50]}...': {e}")

    if not RDKIT_AVAILABLE:
        msg = "RDKit not available. NMR filtering skipped."
        report_progress(msg)
        return candidate_cids, msg # Return original CIDs and a warning

    if required_signals is None and required_ratio is None:
        report_progress("No NMR constraints specified. Skipping NMR filtering.")
        return candidate_cids, None # No filtering needed

    if not candidate_cids:
        report_progress("No candidate CIDs to filter by NMR.")
        return [], None # No candidates to filter

    report_progress(f"Starting NMR filtering for {len(candidate_cids)} CIDs.")
    report_progress(f"  Required Signals: {required_signals if required_signals is not None else 'Any'}")
    report_progress(f"  Required Ratio: {':'.join(map(str,required_ratio)) if required_ratio else 'Any'}")
    
    # Step 1: Fetch SMILES for all candidates
    # The progress_callback is passed down to fetch_smiles_batch
    cid_smiles_map, smiles_fetch_warning = fetch_smiles_batch(candidate_cids, progress_callback)
    
    if not cid_smiles_map:
        msg = "No valid SMILES obtained for any candidate CIDs. Cannot perform NMR filtering."
        report_progress(msg)
        # Concatenate warnings if smiles_fetch_warning exists
        final_warning = msg
        if smiles_fetch_warning:
            final_warning = f"{smiles_fetch_warning}; {msg}"
        return [], final_warning

    # Step 2: Process each compound
    report_progress(f"Predicting NMR properties for {len(cid_smiles_map)} compounds with SMILES...")
    
    filtered_cids: List[int] = []
    nmr_prediction_errors = 0
    total_to_process = len(cid_smiles_map)
    processed_count = 0

    for cid, smiles_string in cid_smiles_map.items():
        processed_count += 1
        if processed_count % 200 == 1 or processed_count == total_to_process : # Log progress periodically
             report_progress(f"  NMR prediction progress: {processed_count}/{total_to_process}...")

        num_signals, integrations_dict = get_hydrogen_environments(smiles_string)

        if num_signals is None: # Error during RDKit processing for this SMILES
            nmr_prediction_errors += 1
            continue # Skip this compound

        # Filter by number of signals
        if required_signals is not None and num_signals != required_signals:
            logger.debug(f"CID {cid}: Predicted signals {num_signals} != required {required_signals}. Filtered out.")
            continue

        # Filter by ratio (if required)
        if required_ratio is not None:
            predicted_ratio = calculate_predicted_ratio(integrations_dict)
            if predicted_ratio is None: # Error in ratio calculation
                logger.debug(f"CID {cid}: Error calculating predicted ratio. Filtered out.")
                nmr_prediction_errors +=1 # Count as an error if ratio was required but couldn't be calculated
                continue
            # Compare sorted lists for ratios
            if predicted_ratio != required_ratio:
                logger.debug(f"CID {cid}: Predicted ratio {predicted_ratio} != required {required_ratio}. Filtered out.")
                continue
        
        # If all checks passed
        filtered_cids.append(cid)

    report_progress(f"NMR filtering complete. {len(filtered_cids)} CIDs matched criteria.")

    # Compile warning messages
    warnings_list = []
    if smiles_fetch_warning:
        warnings_list.append(smiles_fetch_warning)
    if nmr_prediction_errors > 0:
        warnings_list.append(f"{nmr_prediction_errors} compound(s) skipped due to NMR prediction errors.")
    
    final_warning_message = "; ".join(warnings_list) if warnings_list else None
    
    return filtered_cids, final_warning_message