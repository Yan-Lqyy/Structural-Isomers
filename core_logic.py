# core_logic.py

import requests
import json
import time
import math
import logging
from collections import Counter
# *** Ensure Callable is imported ***
from typing import List, Dict, Optional, Any, Tuple, Callable

# --- RDKit Import ---
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

# --- Setup Logging ---
# Match app.py level or set independently. DEBUG is good for development.
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# --- Constants ---
BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
REQUEST_TIMEOUT = 120 # seconds
SMILES_BATCH_SIZE = 500

# --- SMARTS Library Data ---
# (Keep the full SMARTS_LIBRARY_DATA dictionary here)
SMARTS_LIBRARY_DATA = {
    # Carbon Types
    'alkyl_carbon': {'name': 'Alkyl Carbon (sp3)', 'smarts': '[CX4]', 'category': 'Carbon Types', 'desc': 'Saturated carbon.'},
    'allenic_carbon': {'name': 'Allenic Carbon', 'smarts': '[$([CX2](=C)=C)]', 'category': 'Carbon Types', 'desc': 'Center carbon in C=C=C.'},
    'vinylic_carbon': {'name': 'Vinylic Carbon (sp2 alkene)', 'smarts': '[$([CX3]=[CX3])]', 'category': 'Carbon Types', 'desc': 'Carbon in a C=C double bond.'},
    'acetylenic_carbon': {'name': 'Acetylenic Carbon (sp alkyne)', 'smarts': '[$([CX2]#C)]', 'category': 'Carbon Types', 'desc': 'Carbon in a C#C triple bond.'},
    'aromatic_carbon': {'name': 'Aromatic Carbon', 'smarts': 'c', 'category': 'Carbon Types', 'desc': 'Carbon in an aromatic ring.'},
    # Carbonyl and Related (C & O)
    'carbonyl_generic': {'name': 'Carbonyl (Generic, C=O)', 'smarts': '[CX3]=[OX1]', 'category': 'Carbonyl Groups', 'desc': 'Low specificity C=O group.'},
    'carbonyl_reson': {'name': 'Carbonyl (Resonance)', 'smarts': '[$([CX3]=[OX1]),$([CX3+]-[OX1-])]', 'category': 'Carbonyl Groups', 'desc': 'Hits C=O or C+-O- resonance.'},
    'aldehyde': {'name': 'Aldehyde (-CHO)', 'smarts': '[CX3H1](=O)[#6]', 'category': 'Carbonyl Groups', 'desc': 'Requires H on carbonyl C.'},
    'ketone': {'name': 'Ketone (R-CO-R)', 'smarts': '[#6][CX3](=O)[#6]', 'category': 'Carbonyl Groups', 'desc': 'Carbonyl flanked by two carbons.'},
    'anhydride': {'name': 'Anhydride (Acyl)', 'smarts': '[CX3](=[OX1])[OX2][CX3](=[OX1])', 'category': 'Carbonyl Groups', 'desc': 'R-CO-O-CO-R group.'},
    'amide': {'name': 'Amide (R-CONH-R)', 'smarts': '[NX3][CX3](=[OX1])[#6]', 'category': 'Carbonyl Groups', 'desc': 'Amide group (C-version).'},
    'carbamate_ester': {'name': 'Carbamate Ester', 'smarts': '[NX3][CX3](=[OX1])[OX2H0]', 'category': 'Carbonyl Groups', 'desc': 'R-NH-CO-O-R group.'},
    'carbamate_acid': {'name': 'Carbamic Acid/Zwitterion', 'smarts': '[NX3,NX4+][CX3](=[OX1])[OX2H,OX1-]', 'category': 'Carbonyl Groups', 'desc': 'R-NH-COOH or zwitterion.'},
    'carboxylate_ion': {'name': 'Carboxylate Ion (-COO-)', 'smarts': '[CX3](=O)[O-]', 'category': 'Carbonyl Groups', 'desc': 'Anion of carboxylic/carbamic/carbonic.'},
    'carbonic_ester': {'name': 'Carbonic Ester (Diester)', 'smarts': 'C[OX2][CX3](=[OX1])[OX2]C', 'category': 'Carbonyl Groups', 'desc': 'RO-CO-OR group.'},
    'carboxylic_acid': {'name': 'Carboxylic Acid (-COOH)', 'smarts': '[CX3](=O)[OX2H1]', 'category': 'Carbonyl Groups', 'desc': 'COOH group.'},
    'carboxylic_acid_ion': {'name': 'Carboxylic Acid or Ion', 'smarts': '[CX3](=O)[OX1H0-,OX2H1]', 'category': 'Carbonyl Groups', 'desc': 'COOH or COO-.'},
    'ester': {'name': 'Ester (R-COO-R)', 'smarts': '[#6][CX3](=O)[OX2H0][#6]', 'category': 'Carbonyl Groups', 'desc': 'Ester group; won\'t hit formic anhydride.'},
     # Ether
    'ether': {'name': 'Ether (R-O-R)', 'smarts': '[OD2]([#6])[#6]', 'category': 'Oxygen Groups', 'desc': 'Oxygen bonded to two carbons.'},
    # Hydroxyl
    'hydroxyl': {'name': 'Hydroxyl (-OH)', 'smarts': '[OX2H]', 'category': 'Oxygen Groups', 'desc': 'Generic -OH group.'},
    'alcohol': {'name': 'Alcohol (Aliphatic OH)', 'smarts': '[#6][OX2H]', 'category': 'Oxygen Groups', 'desc': '-OH on sp3 carbon.'},
    'enol': {'name': 'Enol (C=C-OH)', 'smarts': '[OX2H][#6X3]=[#6]', 'category': 'Oxygen Groups', 'desc': '-OH on sp2 alkene carbon.'},
    'phenol': {'name': 'Phenol (Ar-OH)', 'smarts': '[OX2H][cX3]:[c]', 'category': 'Oxygen Groups', 'desc': '-OH on aromatic carbon.'},
    'hydroxyl_acidic': {'name': 'Acidic Hydroxyl', 'smarts': '[$([OH]-*=[!#6])]', 'category': 'Oxygen Groups', 'desc': 'OH bonded to atom multiple-bonded to heteroatom.'},
    # Peroxide
    'peroxide': {'name': 'Peroxide (R-O-O-R)', 'smarts': '[OX2,OX1-][OX2,OX1-]', 'category': 'Oxygen Groups', 'desc': 'O-O bond, includes anions.'},
    # Nitrogen Groups (excluding Amide, Carbamate - listed under Carbonyl)
    'amine_primary': {'name': 'Primary Amine (R-NH2)', 'smarts': '[NX3;H2;!$(NC=[!#6]);!$(NC#[!#6])][#6]', 'category': 'Nitrogen Groups', 'desc': 'Primary amine, not amide/cyanamide.'},
    'amine_secondary': {'name': 'Secondary Amine (R-NH-R)', 'smarts': '[NX3;H1;!$(NC=O)][#6]', 'category': 'Nitrogen Groups', 'desc': 'Secondary amine, not amide.'},
    'amine_tertiary': {'name': 'Tertiary Amine (R3N)', 'smarts': '[NX3;H0;!$(NC=O)]([#6])[#6]', 'category': 'Nitrogen Groups', 'desc': 'Tertiary amine, not amide.'},
    'amine_prim_sec': {'name': 'Primary or Secondary Amine', 'smarts': '[NX3;H2,H1;!$(NC=O)]', 'category': 'Nitrogen Groups', 'desc': 'Primary or secondary, not amide.'},
    'enamine': {'name': 'Enamine (C=C-N)', 'smarts': '[NX3][CX3]=[CX3]', 'category': 'Nitrogen Groups', 'desc': 'Nitrogen attached to C=C.'},
    'aniline_enamine_N': {'name': 'Aniline or Enamine Nitrogen', 'smarts': '[NX3][$(C=C),$(cc)]', 'category': 'Nitrogen Groups', 'desc': 'Nitrogen attached to C=C or aromatic ring.'},
    'azide_group': {'name': 'Azide Group (-N3)', 'smarts': '[$(*-[NX2-]-[NX2+]#[NX1]),$(*-[NX2]=[NX2+]=[NX1-])]', 'category': 'Nitrogen Groups', 'desc': 'Attached azide group.'},
    'azide_ion': {'name': 'Azide Ion (N3-)', 'smarts': '[$([NX1-]=[NX2+]=[NX1-]),$([NX1]#[NX2+]-[NX1-2])]', 'category': 'Nitrogen Groups', 'desc': 'Azide anion.'},
    'azo_diazene': {'name': 'Azo / Diazene (N=N)', 'smarts': '[NX2]=[NX2]', 'category': 'Nitrogen Groups', 'desc': 'R-N=N-R group.'},
    'azoxy': {'name': 'Azoxy (R-N=N+(O-)-R)', 'smarts': '[$([NX2]=[NX3+]([O-])[#6]),$([NX2]=[NX3+0](=[O])[#6])]', 'category': 'Nitrogen Groups', 'desc': 'Azoxy group.'},
    'diazo': {'name': 'Diazo (R=N+=N-)', 'smarts': '[$([#6]=[N+]=[N-]),$([#6-]-[N+]#[N])]', 'category': 'Nitrogen Groups', 'desc': 'Diazo group.'},
    'hydrazine': {'name': 'Hydrazine (R2N-NR2)', 'smarts': '[NX3][NX3]', 'category': 'Nitrogen Groups', 'desc': 'N-N single bond.'},
    'hydrazone': {'name': 'Hydrazone (C=N-NR2)', 'smarts': '[NX3][NX2]=[*]', 'category': 'Nitrogen Groups', 'desc': 'Hydrazone group.'},
    'imine': {'name': 'Imine (C=N)', 'smarts': '[$([CX3]([#6])[#6]),$([CX3H][#6])]=[$([NX2][#6]),$([NX2H])]', 'category': 'Nitrogen Groups', 'desc': 'Substituted or unsubstituted imine.'},
    'iminium': {'name': 'Iminium (C=N+)', 'smarts': '[NX3+]=[CX3]', 'category': 'Nitrogen Groups', 'desc': 'Iminium ion.'},
    'nitrate': {'name': 'Nitrate Group/Ion (-ONO2)', 'smarts': '[$([NX3](=[OX1])(=[OX1])O),$([NX3+]([OX1-])(=[OX1])O)]', 'category': 'Nitrogen Groups', 'desc': 'Nitrate group or anion.'},
    'nitrile': {'name': 'Nitrile (-CN)', 'smarts': '[NX1]#[CX2]', 'category': 'Nitrogen Groups', 'desc': 'Cyano group.'},
    'isonitrile': {'name': 'Isonitrile (-NC)', 'smarts': '[CX1-]#[NX2+]', 'category': 'Nitrogen Groups', 'desc': 'Isocyano group.'},
    'nitro': {'name': 'Nitro (-NO2)', 'smarts': '[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8]', 'category': 'Nitrogen Groups', 'desc': 'Nitro group R-NO2.'},
    'nitroso': {'name': 'Nitroso (-NO)', 'smarts': '[NX2]=[OX1]', 'category': 'Nitrogen Groups', 'desc': 'Nitroso group R-N=O.'},
    'n_oxide': {'name': 'N-Oxide', 'smarts': '[$([#7+][OX1-]),$([#7v5]=[OX1]);!$([#7](~[O])~[O]);!$([#7]=[#7])]', 'category': 'Nitrogen Groups', 'desc': 'N->O coordinate bond, excludes nitro etc.'},
    # Sulfur Groups
    'thiol': {'name': 'Thiol (-SH)', 'smarts': '[#16X2H]', 'category': 'Sulfur Groups', 'desc': 'Sulfhydryl or mercapto group.'},
    'thioamide': {'name': 'Thioamide (R-CSNH-R)', 'smarts': '[NX3][CX3]=[SX1]', 'category': 'Sulfur Groups', 'desc': 'Thioamide group.'},
    'sulfide': {'name': 'Sulfide (R-S-R)', 'smarts': '[#16X2H0][!#16]', 'category': 'Sulfur Groups', 'desc': 'Thioether, C-S-C, excludes disulfides.'},
    'disulfide': {'name': 'Disulfide (R-S-S-R)', 'smarts': '[#16X2H0][#16X2H0]', 'category': 'Sulfur Groups', 'desc': 'S-S bond.'},
    'sulfoxide': {'name': 'Sulfoxide (R-SO-R)', 'smarts': '[$([#16X3](=[OX1])([#6])[#6]),$([#16X3+]([OX1-])([#6])[#6])]', 'category': 'Sulfur Groups', 'desc': 'Carbo-sulfoxide R-S(=O)-R.'},
    'sulfone': {'name': 'Sulfone (R-SO2-R)', 'smarts': '[$([#16X4](=[OX1])(=[OX1])([#6])[#6]),$([#16X4+2]([OX1-])([OX1-])([#6])[#6])]', 'category': 'Sulfur Groups', 'desc': 'Carbo-sulfone R-S(=O)2-R.'},
    'sulfonic_acid': {'name': 'Sulfonic Acid/Ion (-SO3H)', 'smarts': '[$([#16X4](=[OX1])(=[OX1])([#6])[OX2H,OX1H0-]),$([#16X4+2]([OX1-])([OX1-])([#6])[OX2H,OX1H0-])]', 'category': 'Sulfur Groups', 'desc': 'R-SO3H or R-SO3-.'},
    'sulfonate_ester': {'name': 'Sulfonate Ester (-SO3R)', 'smarts': '[$([#16X4](=[OX1])(=[OX1])([#6])[OX2H0]),$([#16X4+2]([OX1-])([OX1-])([#6])[OX2H0])]', 'category': 'Sulfur Groups', 'desc': 'R-SO2-OR group.'},
    'sulfonamide': {'name': 'Sulfonamide (-SO2NR2)', 'smarts': '[$([#16X4]([NX3])(=[OX1])(=[OX1])[#6]),$([#16X4+2]([NX3])([OX1-])([OX1-])[#6])]', 'category': 'Sulfur Groups', 'desc': 'Carbo-sulfonamide R-SO2-NR2.'},
    'sulfate_monoester': {'name': 'Sulfate (Monoester/Acid)', 'smarts': '[$([#16X4](=[OX1])(=[OX1])([OX2H,OX1H0-])[OX2][#6]),$([#16X4+2]([OX1-])([OX1-])([OX2H,OX1H0-])[OX2][#6])]', 'category': 'Sulfur Groups', 'desc': 'RO-SO3H or RO-SO3-.'},
    'sulfate_diester': {'name': 'Sulfate (Diester)', 'smarts': '[$([#16X4](=[OX1])(=[OX1])([OX2][#6])[OX2][#6]),$([#16X4](=[OX1])(=[OX1])([OX2][#6])[OX2][#6])]', 'category': 'Sulfur Groups', 'desc': 'RO-SO2-OR group.'},
    # Halogens
    'halogen': {'name': 'Halogen Atom', 'smarts': '[F,Cl,Br,I]', 'category': 'Halogens', 'desc': 'Any halogen atom.'},
    'alkyl_halide': {'name': 'Alkyl Halide (C-X)', 'smarts': '[#6][F,Cl,Br,I]', 'category': 'Halogens', 'desc': 'Halogen attached to any carbon.'},
    'acyl_halide': {'name': 'Acyl Halide (R-COX)', 'smarts': '[CX3](=[OX1])[F,Cl,Br,I]', 'category': 'Halogens', 'desc': 'Acid halide.'},
     # Phosphorus
    'phosphoric_acid': {'name': 'Phosphoric Acid Group', 'smarts': '[$(P(=[OX1])([$([OX2H]),$([OX1-]),$([OX2]P)])([$([OX2H]),$([OX1-]),$([OX2]P)])[$([OX2H]),$([OX1-]),$([OX2]P)]),$([P+]([OX1-])([$([OX2H]),$([OX1-]),$([OX2]P)])([$([OX2H]),$([OX1-]),$([OX2]P)])[$([OX2H]),$([OX1-]),$([OX2]P)])]', 'category': 'Phosphorus Groups', 'desc': 'Orthophosphoric or polyphosphoric acid/anhydride.'},
    'phosphoric_ester': {'name': 'Phosphoric Ester Group', 'smarts': '[$(P(=[OX1])([OX2][#6])([$([OX2H]),$([OX1-]),$([OX2][#6])])[$([OX2H]),$([OX1-]),$([OX2][#6]),$([OX2]P)]),$([P+]([OX1-])([OX2][#6])([$([OX2H]),$([OX1-]),$([OX2][#6])])[$([OX2H]),$([OX1-]),$([OX2][#6]),$([OX2]P)])]', 'category': 'Phosphorus Groups', 'desc': 'Phosphate ester P-O-C bond.'},
}

# --- Helper to get categories ---
def get_smarts_categories(library_data):
    categories = {}
    for key, data in library_data.items():
        cat = data.get('category', 'Uncategorized')
        if cat not in categories: categories[cat] = []
        categories[cat].append({
            'key': key, 'name': data.get('name', key), 'desc': data.get('desc', '')
        })
    for cat in categories:
        categories[cat].sort(key=lambda x: x['name']) # Sort items within category by name
    # Sort categories themselves by name
    return dict(sorted(categories.items()))

# --- API Call Helper ---
def run_pubchem_search(url: str, params: Optional[Dict] = None, request_type: str = "GET", data: Optional[Dict] = None, retries: int = 2, delay: float = 2.0) -> Optional[Dict]:
    """Makes request to PubChem PUG REST API with retries and error handling."""
    headers = {'Accept': 'application/json', 'User-Agent': 'IsomerWebApp/1.2'} # Update version?
    last_exception = None
    last_status_code = None
    logger.debug(f"PubChem Request ({request_type}): {url} | Params: {params}")
    for attempt in range(retries + 1):
        try:
            if request_type.upper() == "GET":
                response = requests.get(url, params=params, headers=headers, timeout=REQUEST_TIMEOUT)
            elif request_type.upper() == "POST":
                 headers['Content-Type'] = 'application/x-www-form-urlencoded'
                 response = requests.post(url, data=data, headers=headers, timeout=REQUEST_TIMEOUT)
            else:
                 logger.error(f"Unsupported request type '{request_type}'")
                 return {'_error_status': 400, '_error_message': f"Internal Error: Unsupported request type '{request_type}'"}

            last_status_code = response.status_code
            response.raise_for_status() # Raises HTTPError for 4xx/5xx responses
            result_data = response.json()
            logger.debug(f"PubChem request successful (Status: {last_status_code})")
            return result_data

        except requests.exceptions.Timeout as e:
            last_exception = e
            logger.warning(f"PubChem request timed out (Attempt {attempt+1}/{retries+1}). URL: {url}")
        except requests.exceptions.HTTPError as e:
            last_exception = e
            status_code = response.status_code
            server_error_message = f"PubChem API Error {status_code}"
            try: # Try getting detailed error from PubChem JSON response
                error_details = response.json()
                fault_msg = error_details.get('Fault', {}).get('Message')
                details = error_details.get('Fault', {}).get('Details')
                if fault_msg: server_error_message += f": {fault_msg}"
                if details: server_error_message += f" ({'; '.join(details)})"
            except json.JSONDecodeError: server_error_message += f". Response: {response.text[:200]}" # Show partial text if not JSON
            except Exception: pass # Ignore errors getting more details
            if 500 <= status_code <= 599 and attempt < retries: # Retry on 5xx server errors
                 logger.warning(f"PubChem API returned HTTP {status_code} (Attempt {attempt+1}/{retries+1}). Retrying... URL: {url}")
            else: # Non-retryable error (4xx) or last attempt failed
                 logger.error(f"PubChem API returned non-retryable HTTP Error {status_code} for URL {url} (Final Attempt). Message: {server_error_message}")
                 return {'_error_status': status_code, '_error_message': server_error_message}
        except json.JSONDecodeError as e:
            last_exception = e
            logger.error(f"Failed to decode JSON response from PubChem for {url} (Attempt {attempt+1}/{retries+1}). Status: {last_status_code}. Response: {response.text[:200] if response else 'N/A'}")
            # Treat as internal error as we expect JSON
            return {'_error_status': 500, '_error_message': "Invalid JSON response from PubChem"}
        except requests.exceptions.RequestException as e:
             last_exception = e
             logger.error(f"Network or request error occurred: {e} (Attempt {attempt+1}/{retries+1}). URL: {url}")
        except Exception as e: # Catch any other unexpected errors
             last_exception = e
             logger.exception(f"An unexpected error occurred during API call to {url} (Attempt {attempt+1}/{retries+1}).") # Log full trace for unexpected

        # If we are here, an error occurred and we might retry
        if attempt < retries:
            calculated_delay = delay * (1.5**attempt) # Exponential backoff
            logger.info(f"Waiting {calculated_delay:.2f}s before retry...")
            time.sleep(calculated_delay)
        else:
            break # Max retries reached

    # Fallback after all retries failed
    error_message = f"Request failed after {retries+1} attempts. Last status: {last_status_code}. Last error: {type(last_exception).__name__ if last_exception else 'N/A'}"
    logger.error(f"{error_message} URL: {url}")
    return {'_error_status': last_status_code or 500, '_error_message': error_message}


# --- Structure Search Function (Depends only on run_pubchem_search) ---
def find_candidates_by_structure(
    formula: str,
    constraints_dict: Dict[str, str],
    # *** ADD progress_callback argument ***
    progress_callback: Optional[Callable[[str], None]] = None
) -> Tuple[Optional[List[int]], Optional[str]]:
    """
    Finds candidate CIDs by formula and structural constraints, reporting progress.
    """
    def report_progress(message: str):
        """Helper to log and call callback if available."""
        logger.info(message) # Log the message regardless
        if progress_callback:
            try:
                progress_callback(message) # Send to client via callback
            except Exception as e:
                 # Log error in callback but don't stop the search process
                 logger.warning(f"Progress callback failed for message '{message[:50]}...': {e}")

    report_progress(f"Starting structure search for formula: '{formula}'")
    formula_url = f"{BASE_URL}/compound/fastformula/{formula}/cids/JSON"
    formula_params = {'list_return': 'cachekey'}
    formula_data = run_pubchem_search(formula_url, params=formula_params)

    if '_error_status' in formula_data:
        msg = formula_data.get('_error_message', 'Unknown error during formula search')
        report_progress(f"Formula search failed: {msg}") # Report failure via callback too
        return None, msg
    if not formula_data or 'IdentifierList' not in formula_data or 'CacheKey' not in formula_data['IdentifierList'] or formula_data['IdentifierList'].get('Size', 0) == 0:
        report_progress(f"No compounds found matching formula '{formula}'.")
        return [], None # Return empty list, no error message needed

    current_cachekey = formula_data['IdentifierList']['CacheKey']
    initial_size = formula_data['IdentifierList']['Size']
    report_progress(f"Formula matched {initial_size} compounds. (CacheKey: {current_cachekey[:10]}...)")

    if not constraints_dict:
        report_progress("No structural constraints specified. Retrieving all candidates from cache...")
        final_cids_url = f"{BASE_URL}/compound/listkey/{current_cachekey}/cids/JSON"
        final_cids_data = run_pubchem_search(final_cids_url)
        if '_error_status' in final_cids_data:
            msg = final_cids_data.get('_error_message', 'Error retrieving CIDs from cache')
            report_progress(f"Failed to retrieve full list: {msg}")
            return None, msg
        if final_cids_data and 'IdentifierList' in final_cids_data and 'CID' in final_cids_data['IdentifierList']:
            cids = final_cids_data['IdentifierList']['CID']
            report_progress(f"Retrieved {len(cids)} candidates matching formula only.")
            return cids, None
        else:
            msg = "Error parsing candidate CIDs from cachekey."
            report_progress(msg)
            return None, msg

    report_progress(f"Applying {len(constraints_dict)} structural constraint(s)...")
    num_constraints = len(constraints_dict)
    constraint_items = list(constraints_dict.items())

    for i, (key, smarts) in enumerate(constraint_items):
        is_last_constraint = (i == num_constraints - 1)
        report_progress(f"  Applying constraint {i+1}/{num_constraints}: '{key}'...")

        refine_url = f"{BASE_URL}/compound/fastsubstructure/smarts/{smarts}/cids/JSON"
        refine_params = {'cachekey': current_cachekey}
        if not is_last_constraint:
            refine_params['list_return'] = 'cachekey'
            expected_key = 'CacheKey'
        else:
            expected_key = 'CID' # Final step returns CIDs

        # Make the API call
        refine_data = run_pubchem_search(refine_url, params=refine_params)

        # Handle API call errors
        if '_error_status' in refine_data:
            msg = refine_data.get('_error_message', f"Error during refinement for '{key}'")
            report_progress(f"Constraint application failed: {msg}")
            return None, msg
        if not refine_data or 'IdentifierList' not in refine_data:
            msg = f"Unexpected empty response during refinement for '{key}'."
            report_progress(msg)
            return None, msg

        # Process successful response
        identifier_list = refine_data['IdentifierList']
        result = identifier_list.get(expected_key)
        size = identifier_list.get('Size', 0)

        if result is None and size == 0:
            report_progress(f"Constraint '{key}' yielded 0 results. No compounds match all constraints.")
            return [], None # Found zero matching all constraints, return empty list
        if result is None:
            # This case means size > 0 but the expected key (CacheKey/CID) is missing - unexpected
            msg = f"Failed to get expected result ('{expected_key}') for constraint '{key}', though size reported as {size}."
            report_progress(msg)
            return None, msg

        # Update state for next loop or return final results
        if not is_last_constraint:
            current_cachekey = result
            report_progress(f"  Constraint '{key}' matched. {size} compounds remaining. (CacheKey: {current_cachekey[:10]}...)")
        else:
            candidate_cids = result
            report_progress(f"Finished structural filtering. Found {len(candidate_cids)} candidate(s).")
            return candidate_cids, None

    # Should not be reached if constraints exist
    msg = "Constraint application logic error (unexpected end of loop)."
    report_progress(msg)
    return None, msg


# *** >>> MOVE NMR HELPER FUNCTIONS HERE <<< ***
# --- NMR Prediction Helpers ---
def get_hydrogen_environments(smiles_string: str) -> Tuple[Optional[int], Optional[Dict[int, int]]]:
    """Calculates number of signals and integration counts from SMILES using RDKit."""
    if not RDKIT_AVAILABLE or not smiles_string:
        if RDKIT_AVAILABLE and not smiles_string: logger.warning("get_hydrogen_environments called with empty SMILES.")
        return None, None
    try:
        mol = Chem.MolFromSmiles(smiles_string)
        if mol is None: logger.debug(f"RDKit failed to parse SMILES: {smiles_string}"); return None, None
        mol_h = Chem.AddHs(mol)
        if mol_h is None: logger.warning(f"RDKit failed to add Hydrogens to SMILES: {smiles_string}"); return None, None

        has_hydrogens = any(atom.GetAtomicNum() == 1 for atom in mol_h.GetAtoms())
        if not has_hydrogens: logger.debug(f"No hydrogens found in: {smiles_string}"); return 0, {}

        ranks = list(Chem.CanonicalRankAtoms(mol_h, breakTies=False))
        h_counts_by_rank = Counter()
        for atom in mol_h.GetAtoms():
            if atom.GetAtomicNum() == 1:
                try: h_counts_by_rank[ranks[atom.GetIdx()]] += 1
                except IndexError: logger.warning(f"Index error H rank {atom.GetIdx()} in {smiles_string}"); continue

        num_signals = len(h_counts_by_rank)
        integration_counts_dict = dict(h_counts_by_rank)
        if has_hydrogens and num_signals == 0: logger.warning(f"RDKit found 0 H environments for {smiles_string} despite H present.")

        # *** Make NMR log message more CONCISE ***
        logger.debug(f"NMR Prediction for {smiles_string[:30]}... : Found {num_signals} signal(s)")
        # logger.debug(f"SMILES {smiles_string}: Found {num_signals} signals, Integrations: {integration_counts_dict}") # Old verbose log
        return num_signals, integration_counts_dict

    except Exception as e:
        logger.warning(f"RDKit error processing SMILES '{smiles_string}': {e}", exc_info=False)
        return None, None

def gcd_list(numbers: List[int]) -> int:
    """Calculates the greatest common divisor of a list of integers."""
    positive_numbers = [int(n) for n in numbers if n > 0]
    if not positive_numbers: return 1
    if len(positive_numbers) == 1: return positive_numbers[0]
    result = positive_numbers[0]
    for i in range(1, len(positive_numbers)): result = math.gcd(result, positive_numbers[i])
    return result if result > 0 else 1

def calculate_predicted_ratio(integration_counts_dict: Optional[Dict[int, int]]) -> Optional[List[int]]:
    """Calculates the simplified NMR ratio from integration counts."""
    if integration_counts_dict is None: return None
    if not integration_counts_dict: return []
    counts = list(integration_counts_dict.values())
    if not counts: return []
    positive_counts = [c for c in counts if c > 0]
    if not positive_counts: logger.warning(f"Ratio calc found no positive counts in {counts}"); return []
    common_divisor = gcd_list(positive_counts)
    if common_divisor <= 0: logger.error(f"Invalid GCD ({common_divisor}) for counts: {positive_counts}"); return None
    ratios = sorted([count // common_divisor for count in counts if count > 0])
    logger.debug(f"Calculated ratio {ratios} from counts {counts} (GCD: {common_divisor})")
    return ratios

def format_nmr_prediction(num_signals: Optional[int], integration_counts_dict: Optional[Dict[int, int]]) -> str:
    """Formats the predicted NMR data into a user-friendly string for display."""
    if num_signals is None: return "N/A (Prediction Error)"
    if num_signals == 0: return "0 Signals (No Equivalent H Expected)"
    integration_str = "N/A"
    if integration_counts_dict:
        counts_display = sorted(list(integration_counts_dict.values()), reverse=True)
        integration_str = ", ".join([f"{count}H" for count in counts_display])
    ratio_list = calculate_predicted_ratio(integration_counts_dict)
    ratio_str = "N/A"
    if ratio_list is None: ratio_str = "Error"
    elif ratio_list: ratio_str = ":".join(map(str, ratio_list))
    # Return string for display on results page
    return f"{num_signals} Signals; Ratio: {ratio_str} ({integration_str})"
# *** >>> END OF MOVED NMR HELPER FUNCTIONS <<< ***


# --- Batch SMILES Fetching (Depends on run_pubchem_search) ---
def fetch_smiles_batch(
    cids: List[int],
    # *** ADD progress_callback argument ***
    progress_callback: Optional[Callable[[str], None]] = None
) -> Tuple[Dict[int, str], Optional[str]]:
    """
    Fetches Isomeric (or fallback Canonical) SMILES for a list of CIDs in batches, reporting progress.
    """
    def report_progress(message: str):
        """Helper to log and call callback if available."""
        logger.info(message)
        if progress_callback:
            try: progress_callback(message)
            except Exception as e: logger.warning(f"Progress callback failed: {e}")

    if not cids: return {}, None

    report_progress(f"Fetching SMILES for {len(cids)} CIDs (Batch size: {SMILES_BATCH_SIZE})...")
    cid_smiles_map = {}
    warnings = []
    fetch_errors = 0
    total_batches = (len(cids) + SMILES_BATCH_SIZE - 1) // SMILES_BATCH_SIZE

    for i in range(0, len(cids), SMILES_BATCH_SIZE):
        batch_num = i // SMILES_BATCH_SIZE + 1
        batch_cids = cids[i:i+SMILES_BATCH_SIZE]
        cids_str = ','.join(map(str, batch_cids))
        properties = "IsomericSMILES,CanonicalSMILES"
        url = f"{BASE_URL}/compound/cid/{cids_str}/property/{properties}/JSON"
        # Report progress less verbosely for many batches
        if batch_num % 5 == 1 or batch_num == total_batches or total_batches <= 5 :
             report_progress(f"  Fetching SMILES batch {batch_num}/{total_batches}...")

        data = run_pubchem_search(url)

        if data is None or '_error_status' in data:
             error_msg = data.get('_error_message', f'API error fetching SMILES batch {batch_num}.') if isinstance(data, dict) else 'Unknown API error'
             logger.error(error_msg) # Log error internally
             # Don't report API errors via progress callback, keep log cleaner
             warnings.append(f"Failed to fetch SMILES batch {batch_num}")
             fetch_errors += len(batch_cids)
             continue

        if 'PropertyTable' in data and 'Properties' in data['PropertyTable']:
            props_found_in_batch = 0
            for prop_entry in data['PropertyTable']['Properties']:
                cid = prop_entry.get('CID')
                if cid not in batch_cids: continue # Safety check
                iso_smiles = prop_entry.get('IsomericSMILES')
                can_smiles = prop_entry.get('CanonicalSMILES')
                valid_smiles = iso_smiles if iso_smiles and iso_smiles.lower() != 'n/a' else can_smiles if can_smiles and can_smiles.lower() != 'n/a' else None
                if valid_smiles: cid_smiles_map[cid] = valid_smiles; props_found_in_batch += 1
            if props_found_in_batch < len(batch_cids):
                 missing_count = len(batch_cids) - props_found_in_batch
                 logger.warning(f"Missing SMILES properties for {missing_count} CIDs in batch {batch_num}.")
                 # Could add to warnings list if needed on client side
        else:
            logger.warning(f"Unexpected response structure for SMILES batch {batch_num}. Data: {str(data)[:200]}")
            warnings.append(f"Unexpected response structure for SMILES batch {batch_num}")
            fetch_errors += len(batch_cids)

    num_missing = len(cids) - len(cid_smiles_map)
    if num_missing > 0:
        miss_msg = f"Could not retrieve valid SMILES for {num_missing} out of {len(cids)} CIDs."
        logger.warning(miss_msg)
        warnings.append(miss_msg)

    warning_message = "; ".join(warnings) if warnings else None
    report_progress(f"Finished fetching SMILES. Retrieved for {len(cid_smiles_map)}/{len(cids)} CIDs.")
    return cid_smiles_map, warning_message


# --- NMR Filtering Logic (Depends on fetch_smiles_batch and NMR Helpers) ---
# *** THIS FUNCTION MUST COME *AFTER* THE NMR HELPER DEFINITIONS ***
def filter_candidates_by_nmr(
    candidate_cids: List[int],
    required_signals: Optional[int],
    required_ratio: Optional[List[int]],
    # *** ADD progress_callback argument ***
    progress_callback: Optional[Callable[[str], None]] = None
) -> Tuple[List[int], Optional[str]]:
    """
    Filters candidates by NMR using batch SMILES fetching, reporting progress.
    """
    def report_progress(message: str):
        """Helper to log and call callback if available."""
        logger.info(message)
        if progress_callback:
            try: progress_callback(message)
            except Exception as e: logger.warning(f"Progress callback failed: {e}")

    if not RDKIT_AVAILABLE: return candidate_cids, "Skipping NMR filtering (RDKit not available)."
    if required_signals is None and required_ratio is None: return candidate_cids, None
    if not candidate_cids: return [], None

    report_progress(f"Starting NMR filtering for {len(candidate_cids)} candidates...")
    report_progress(f"(Required Signals: {required_signals if required_signals is not None else 'Any'}, Required Ratio: {required_ratio if required_ratio else 'Any'})")

    # Step 3a: Fetch SMILES (will report its own progress via callback)
    start_smiles_time = time.time()
    cid_smiles_map, smiles_warning = fetch_smiles_batch(candidate_cids, progress_callback=progress_callback)
    end_smiles_time = time.time()
    report_progress(f"SMILES fetching step took {end_smiles_time - start_smiles_time:.2f}s.")

    if not cid_smiles_map:
        msg = "Cannot perform NMR filtering: No valid SMILES obtained."
        report_progress(msg)
        # Pass smiles_warning back if it exists
        return [], smiles_warning

    # Step 3b: Filter by NMR properties
    report_progress(f"Processing {len(cid_smiles_map)} compounds with SMILES for NMR properties...")
    final_cids = []
    processed_count = 0
    skipped_nmr = 0
    start_filter_time = time.time()
    total_to_filter = len(cid_smiles_map)

    for cid, smiles in cid_smiles_map.items():
        processed_count += 1
        # Log progress less frequently for performance if many compounds
        if processed_count % 500 == 1 or processed_count == total_to_filter:
             elapsed = time.time() - start_filter_time
             report_progress(f"  NMR Filter Progress: {processed_count}/{total_to_filter} compounds processed ({elapsed:.1f}s)...")

        # --- Get NMR prediction ---
        # get_hydrogen_environments logs concise debug messages internally now
        num_signals, integrations_dict = get_hydrogen_environments(smiles)

        if num_signals is None: # Indicates prediction error
            skipped_nmr += 1
            continue # Skip this compound

        # --- Apply Filters ---
        # 1. Signal Count Filter
        if required_signals is not None and num_signals != required_signals:
            continue # Doesn't match signal count requirement

        # 2. Ratio Filter
        if required_ratio is not None:
            predicted_ratio = calculate_predicted_ratio(integrations_dict)
            # Ratio prediction could fail (None), or might not match
            if predicted_ratio is None or predicted_ratio != required_ratio:
                continue # Doesn't match ratio requirement or prediction failed

        # If we reach here, the compound passed all NMR checks
        final_cids.append(cid)

    end_filter_time = time.time()
    report_progress(f"NMR property calculation & filtering took {end_filter_time - start_filter_time:.2f}s.")

    # --- Compile Warnings ---
    filter_warnings = []
    if smiles_warning: filter_warnings.append(smiles_warning)
    if skipped_nmr > 0:
        skip_msg = f"Skipped {skipped_nmr} compounds due to NMR prediction errors during filtering."
        report_progress(skip_msg) # Report this summary via callback
        filter_warnings.append(skip_msg)
    final_warning_msg = "; ".join(filter_warnings) if filter_warnings else None

    report_progress(f"Finished NMR filtering. {len(final_cids)} compounds matched criteria.")
    return final_cids, final_warning_msg


# --- Detail Fetching (Depends on run_pubchem_search and format_nmr_prediction) ---
# *** THIS FUNCTION MUST COME *AFTER* format_nmr_prediction DEFINITION ***
# (No progress callback needed here as it's relatively quick and done after main search)
def fetch_compound_details(
    cids: List[int],
    detail_image_width: int,
    detail_image_height: int
    ) -> Dict[int, Dict[str, Any]]:
    """Fetches details (properties, description, links, NMR prediction) for a list of CIDs."""
    if not cids: return {}
    logger.info(f"Fetching details for top {len(cids)} compound(s)...")
    compound_data = {cid: {} for cid in cids} # Pre-initialize

    # --- Batch Fetch Properties ---
    properties_to_fetch = ["Title", "IUPACName", "IsomericSMILES", "InChI", "InChIKey"]
    cids_str = ','.join(map(str, cids))
    props_url = f"{BASE_URL}/compound/cid/{cids_str}/property/{','.join(properties_to_fetch)}/JSON"
    props_data = run_pubchem_search(props_url)

    if props_data and 'PropertyTable' in props_data and 'Properties' in props_data['PropertyTable']:
        props_found = 0
        for prop_entry in props_data['PropertyTable']['Properties']:
            cid = prop_entry.get('CID')
            if cid in compound_data:
                compound_data[cid].update(prop_entry)
                props_found += 1
        logger.info(f"Fetched properties for {props_found}/{len(cids)} CIDs.")
    elif '_error_status' in props_data:
        logger.warning(f"Could not fetch properties for details: {props_data.get('_error_message')}")
    else:
        logger.warning("Could not parse properties response for details.")

    # --- Batch Fetch Descriptions ---
    desc_url = f"{BASE_URL}/compound/cid/{cids_str}/description/JSON"
    desc_data = run_pubchem_search(desc_url)
    if desc_data and 'InformationList' in desc_data and 'Information' in desc_data['InformationList']:
        desc_found = 0
        # Find the first valid description per CID (often multiple sources exist)
        processed_cids_desc = set()
        for info in desc_data['InformationList']['Information']:
            cid = info.get('CID')
            description = info.get('Description')
            if cid in compound_data and description and cid not in processed_cids_desc:
                compound_data[cid]['Description'] = description
                processed_cids_desc.add(cid)
                desc_found += 1
        logger.info(f"Fetched descriptions for {desc_found}/{len(cids)} CIDs.")
    elif '_error_status' in desc_data:
        logger.warning(f"Could not fetch descriptions for details: {desc_data.get('_error_message')}")
    else:
        logger.warning("Could not parse descriptions response for details.")


    # --- Add Links and Predict NMR ---
    logger.info("Adding links and predicting NMR for detailed display...")
    nmr_predictions_done = 0
    for cid in cids:
        if cid in compound_data:
            # Add links
            compound_data[cid]['PubChemLink'] = f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}"
            compound_data[cid]['ImageLink'] = f"https://pubchem.ncbi.nlm.nih.gov/image/imagefly.cgi?cid={cid}&width={detail_image_width}&height={detail_image_height}"

            # Predict NMR (using SMILES fetched in properties)
            smiles = compound_data[cid].get('IsomericSMILES')
            if RDKIT_AVAILABLE and smiles and smiles.lower() != 'n/a':
                 signals, integrations = get_hydrogen_environments(smiles)
                 compound_data[cid]['PredictedNMR'] = format_nmr_prediction(signals, integrations)
                 if signals is not None: nmr_predictions_done += 1
            elif not RDKIT_AVAILABLE:
                 compound_data[cid]['PredictedNMR'] = "N/A (RDKit Missing)"
            else:
                 compound_data[cid]['PredictedNMR'] = "N/A (No SMILES)"
        else:
            # This shouldn't happen if compound_data was initialized correctly
             logger.warning(f"CID {cid} missing from compound_data during detail processing.")

    logger.info(f"Predicted NMR for {nmr_predictions_done}/{len(cids)} compounds in details.")
    logger.info("Detail fetching complete.")
    return compound_data