# core_logic.py

import requests
import json
import time
import math
import logging
from collections import Counter
from typing import List, Dict, Optional, Any, Tuple

# --- RDKit Import ---
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

# --- Setup Logging ---
# Use basicConfig for simple module logging if not run as main app
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# --- Constants ---
BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
# Timeout for individual PubChem requests
REQUEST_TIMEOUT = 120 # seconds

# --- SMARTS Library Data ---
# (Keep the full SMARTS_LIBRARY_DATA dictionary here as provided in the previous step)
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
    """Groups SMARTS library data by category for display."""
    categories = {}
    for key, data in library_data.items():
        cat = data.get('category', 'Uncategorized')
        if cat not in categories: categories[cat] = []
        categories[cat].append({
            'key': key,
            'name': data.get('name', key),
            'desc': data.get('desc', '')
        })
    # Sort categories alphabetically, then items within categories by key
    sorted_categories = {cat: sorted(items, key=lambda x: x['key'])
                         for cat, items in categories.items()}
    return dict(sorted(sorted_categories.items()))

# --- API Call Helper (with retries) ---
def run_pubchem_search(url: str, params: Optional[Dict] = None, request_type: str = "GET", data: Optional[Dict] = None, retries: int = 2, delay: float = 2.0) -> Optional[Dict]:
    """Runs a generic PubChem PUG REST request with retries for 5xx errors."""
    headers = {'Accept': 'application/json', 'User-Agent': 'IsomerWebApp/1.0'} # Added User-Agent
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
            response.raise_for_status() # Check for HTTP errors

            result_data = response.json()
            logger.debug(f"PubChem request successful (Status: {last_status_code})")
            return result_data # Success

        except requests.exceptions.Timeout as e:
            last_exception = e
            logger.warning(f"PubChem request timed out (Attempt {attempt+1}/{retries+1}). URL: {url}")
            if attempt < retries: time.sleep(delay * (1.5**attempt))
            else: break
        except requests.exceptions.HTTPError as e:
            last_exception = e
            if response is None: # Defensive check
                 logger.error(f"HTTPError but response object is None (Attempt {attempt+1}/{retries+1})")
                 if attempt < retries: time.sleep(delay * (1.5**attempt)); continue
                 else: break

            status_code = response.status_code
            # Extract server message if possible
            server_error_message = f"PubChem API Error {status_code}"
            try:
                error_details = response.json()
                fault_msg = error_details.get('Fault', {}).get('Message')
                details = error_details.get('Fault', {}).get('Details')
                if fault_msg: server_error_message += f": {fault_msg}"
                if details: server_error_message += f" ({'; '.join(details)})"
            except json.JSONDecodeError:
                 server_error_message += f". Response: {response.text[:200]}" # Show start of non-JSON response

            if 500 <= status_code <= 599 and attempt < retries:
                 logger.warning(f"PubChem API returned HTTP {status_code} (Attempt {attempt+1}/{retries+1}). Retrying... URL: {url}")
                 time.sleep(delay * (1.5**attempt))
            else:
                 logger.error(f"PubChem API returned non-retryable HTTP Error {status_code} for URL {url} (Final Attempt). Message: {server_error_message}")
                 return {'_error_status': status_code, '_error_message': server_error_message}
        except json.JSONDecodeError as e:
            last_exception = e
            if response and (not response.text or response.text.isspace()):
                 logger.warning(f"PubChem returned empty/whitespace response (Status: {last_status_code}). URL: {url}")
                 return {} # OK, just no data
            else:
                 logger.error(f"Failed to decode JSON response from PubChem for {url} (Attempt {attempt+1}/{retries+1}). Status: {last_status_code}. Response: {response.text[:200] if response else 'N/A'}")
                 return {'_error_status': 500, '_error_message': "Invalid JSON response from PubChem"}
        except requests.exceptions.RequestException as e:
             last_exception = e
             logger.error(f"Network or request error occurred: {e} (Attempt {attempt+1}/{retries+1}). URL: {url}")
             if attempt < retries: time.sleep(delay * (1.5**attempt))
             else: break
        except Exception as e:
             last_exception = e
             logger.exception(f"An unexpected error occurred during API call to {url} (Attempt {attempt+1}/{retries+1}).") # Use logger.exception
             if attempt < retries: time.sleep(delay * (1.5**attempt))
             else: break

    # If loop finishes without returning, all retries failed
    error_message = f"Request failed after {retries+1} attempts. Last status: {last_status_code}. Last error: {type(last_exception).__name__}"
    logger.error(f"{error_message} URL: {url}")
    return {'_error_status': last_status_code or 500, '_error_message': error_message}


# --- Core Logic: Structure Search ---
def find_candidates_by_structure(formula: str, constraints_dict: Dict[str, str]) -> Tuple[Optional[List[int]], Optional[str]]:
    """Finds candidate CIDs matching formula and structural constraints. Returns (cids, error_message)."""
    logger.info(f"Step 1: Searching PubChem for formula '{formula}'...")
    formula_url = f"{BASE_URL}/compound/fastformula/{formula}/cids/JSON"
    formula_params = {'list_return': 'cachekey'}
    formula_data = run_pubchem_search(formula_url, params=formula_params)

    if '_error_status' in formula_data:
        error_msg = formula_data.get('_error_message', 'Failed to connect to PubChem for formula search.')
        logger.error(f"Formula search failed: {error_msg}")
        return None, error_msg
    if not formula_data or 'IdentifierList' not in formula_data or 'CacheKey' not in formula_data['IdentifierList'] or formula_data['IdentifierList'].get('Size', 0) == 0:
        logger.info(f"Found 0 compounds matching formula '{formula}' or failed to get initial cachekey.")
        return [], None # No candidates found

    current_cachekey = formula_data['IdentifierList']['CacheKey']
    initial_size = formula_data['IdentifierList']['Size']
    logger.info(f"Found initial set of {initial_size} compounds. CacheKey: {current_cachekey[:10]}...")

    # Apply Structural Constraints
    if not constraints_dict:
        logger.info("Step 2: No structural constraints. Retrieving all candidates from initial cachekey...")
        final_cids_url = f"{BASE_URL}/compound/listkey/{current_cachekey}/cids/JSON" # Use listkey endpoint to retrieve from cachekey
        final_cids_data = run_pubchem_search(final_cids_url) # No params needed

        if '_error_status' in final_cids_data:
            error_msg = final_cids_data.get('_error_message', 'Error retrieving candidate CIDs from cachekey.')
            logger.error(error_msg)
            return None, error_msg
        elif final_cids_data and 'IdentifierList' in final_cids_data and 'CID' in final_cids_data['IdentifierList']:
            cids = final_cids_data['IdentifierList']['CID']
            logger.info(f"Retrieved {len(cids)} candidate CIDs.")
            return cids, None
        else:
            logger.error(f"Unexpected structure retrieving CIDs from listkey {current_cachekey[:10]}. Data: {str(final_cids_data)[:200]}")
            return None, "Error parsing candidate CIDs from cachekey."

    logger.info("Step 2: Applying structural constraints...")
    num_constraints = len(constraints_dict)
    constraint_items = list(constraints_dict.items())

    for i, (key, smarts) in enumerate(constraint_items):
        is_last_constraint = (i == num_constraints - 1)
        logger.info(f"Applying constraint {i+1}/{num_constraints}: '{key}'...")

        refine_url = f"{BASE_URL}/compound/fastsubstructure/smarts/{smarts}/cids/JSON"
        refine_params = {'cachekey': current_cachekey}
        if not is_last_constraint:
            refine_params['list_return'] = 'cachekey'
            expected_key = 'CacheKey'
        else:
            expected_key = 'CID'

        refine_data = run_pubchem_search(refine_url, params=refine_params)

        if '_error_status' in refine_data:
             error_msg = refine_data.get('_error_message', f'Error during refinement for constraint {key}.')
             logger.error(error_msg)
             return None, error_msg
        if not refine_data or 'IdentifierList' not in refine_data:
             error_msg = f"Unexpected response structure during refinement for '{key}'. Data: {str(refine_data)[:200]}"
             logger.error(error_msg)
             return None, error_msg

        identifier_list = refine_data['IdentifierList']
        result = identifier_list.get(expected_key)
        size = identifier_list.get('Size', 0)

        if result is None and size == 0:
             logger.info(f"Constraint '{key}' resulted in 0 matches. No candidates remain.")
             return [], None # No candidates left

        if result is None:
            error_msg = f"Failed to get expected result ({expected_key}) for constraint '{key}'. Data: {identifier_list}"
            logger.error(error_msg)
            return None, error_msg

        if not is_last_constraint:
            current_cachekey = result
            logger.info(f"  {size} compounds remain. New CacheKey: {current_cachekey[:10]}...")
        else:
            candidate_cids = result
            logger.info(f"Found {len(candidate_cids)} candidates matching all structural constraints.")
            return candidate_cids, None

    # Fallback
    return None, "Constraint application logic error."


# --- NMR Prediction Functions ---
def get_smiles_from_cid(cid: int) -> Optional[str]:
    """Fetches the Isomeric or Canonical SMILES string."""
    # Use the simpler run_pubchem_search but expect TXT, not JSON
    # This needs adjustment or a separate simple request function
    url_iso = f"{BASE_URL}/compound/cid/{cid}/property/IsomericSMILES/TXT"
    url_canon = f"{BASE_URL}/compound/cid/{cid}/property/CanonicalSMILES/TXT"
    smiles = None
    try:
        response = requests.get(url_iso, timeout=REQUEST_TIMEOUT/4, headers={'User-Agent': 'IsomerWebApp/1.0'})
        response.raise_for_status()
        smiles = response.text.strip()
        if not smiles or smiles.lower() == 'n/a':
            smiles = None # Reset if invalid
    except requests.exceptions.RequestException:
        pass # Try canonical if isomeric fails

    if smiles is None:
        try:
            response = requests.get(url_canon, timeout=REQUEST_TIMEOUT/4, headers={'User-Agent': 'IsomerWebApp/1.0'})
            response.raise_for_status()
            smiles = response.text.strip()
            if not smiles or smiles.lower() == 'n/a':
                 logger.warning(f"Could not retrieve valid SMILES (Isomeric or Canonical) for CID {cid}.")
                 return None
        except requests.exceptions.RequestException:
            logger.warning(f"Failed to retrieve any SMILES for CID {cid}.")
            return None
    return smiles

def get_hydrogen_environments(smiles_string: str) -> Tuple[Optional[int], Optional[Dict[int, int]]]:
    """Determines unique H environments using RDKit."""
    if not RDKIT_AVAILABLE or not smiles_string: return None, None
    try:
        mol = Chem.MolFromSmiles(smiles_string)
        if mol is None: return None, None
        mol_h = Chem.AddHs(mol)
        if mol_h is None: return None, None
        has_hydrogens = any(atom.GetAtomicNum() == 1 for atom in mol_h.GetAtoms())
        if not has_hydrogens: return 0, {}
        ranks = list(Chem.CanonicalRankAtoms(mol_h, breakTies=False))
        h_counts_by_rank = Counter()
        for atom in mol_h.GetAtoms():
            if atom.GetAtomicNum() == 1:
                try: h_counts_by_rank[ranks[atom.GetIdx()]] += 1
                except IndexError: continue # Should be rare
        num_signals = len(h_counts_by_rank)
        integration_counts_dict = dict(h_counts_by_rank)
        if has_hydrogens and num_signals == 0 and not integration_counts_dict:
            logger.warning(f"RDKit found 0 H environments for {smiles_string} despite H atoms present.")
            return 0, {}
        return num_signals, integration_counts_dict
    except Exception as e:
        logger.warning(f"RDKit error processing SMILES {smiles_string}: {e}")
        return None, None

def gcd_list(numbers: List[int]) -> int:
    positive_numbers = [int(n) for n in numbers if n > 0]
    if not positive_numbers: return 1
    result = positive_numbers[0]
    for i in range(1, len(positive_numbers)): result = math.gcd(result, positive_numbers[i])
    return result

def calculate_predicted_ratio(integration_counts_dict: Optional[Dict[int, int]]) -> Optional[List[int]]:
    if not integration_counts_dict: return []
    counts = list(integration_counts_dict.values())
    if not counts: return []
    common_divisor = gcd_list(counts)
    if common_divisor <= 0: return None
    ratios = sorted([count // common_divisor for count in counts])
    return ratios

def format_nmr_prediction(num_signals: Optional[int], integration_counts_dict: Optional[Dict[int, int]]) -> str:
    """Formats NMR prediction string, used by fetch_compound_details."""
    if num_signals is None or integration_counts_dict is None: return "N/A (Error)"
    if num_signals == 0: return "0 Signals (No H)"
    counts = sorted(list(integration_counts_dict.values()), reverse=True)
    integration_str = ", ".join([f"{count}H" for count in counts])
    ratio_list = calculate_predicted_ratio(integration_counts_dict)
    ratio_str = ":".join(map(str, ratio_list)) if ratio_list is not None else "Error"
    # Shorter format for web display within details
    return f"{num_signals} Signals; Ratio: {ratio_str} ({integration_str})"


# --- NMR Filtering Logic ---
def filter_candidates_by_nmr(
    candidate_cids: List[int],
    required_signals: Optional[int],
    required_ratio: Optional[List[int]]
) -> Tuple[List[int], Optional[str]]:
    """Filters candidates by NMR. Returns (filtered_cids, warning_message)."""
    if not RDKIT_AVAILABLE:
        return candidate_cids, "Skipping NMR filtering (RDKit not available)."
    if required_signals is None and required_ratio is None:
        return candidate_cids, None # No NMR constraints

    logger.info(f"Step 3: Filtering {len(candidate_cids)} candidates by NMR constraints...")
    final_cids = []
    processed_count = 0
    skipped_smiles = 0
    skipped_nmr = 0
    start_time = time.time()

    for cid in candidate_cids:
        processed_count += 1
        if processed_count % 200 == 0: # Log progress less often
             elapsed = time.time() - start_time
             logger.info(f"  NMR Filter Progress: {processed_count}/{len(candidate_cids)} ({elapsed:.1f}s)...")

        smiles = get_smiles_from_cid(cid)
        if not smiles: skipped_smiles += 1; continue

        num_signals, integrations_dict = get_hydrogen_environments(smiles)
        if num_signals is None or integrations_dict is None: skipped_nmr += 1; continue

        # Apply checks
        signal_match = (required_signals is None) or (num_signals == required_signals)
        if not signal_match: continue

        ratio_match = True # Assume match if no ratio required
        if required_ratio is not None:
            predicted_ratio = calculate_predicted_ratio(integrations_dict)
            ratio_match = (predicted_ratio is not None and predicted_ratio == required_ratio)
        if not ratio_match: continue

        # Passed all checks
        final_cids.append(cid)

    end_time = time.time()
    logger.info(f"NMR Filtering Complete in {end_time - start_time:.2f}s.")
    warning_msg = None
    warnings = []
    if skipped_smiles > 0: warnings.append(f"skipped {skipped_smiles} due to SMILES errors")
    if skipped_nmr > 0: warnings.append(f"skipped {skipped_nmr} due to NMR prediction errors")
    if warnings: warning_msg = f"NMR Filtering Warnings: {'; '.join(warnings)}."

    logger.info(f"Found {len(final_cids)} compounds matching all criteria after NMR filter.")
    return final_cids, warning_msg

# --- Detail Fetching ---
def fetch_compound_details(
    cids: List[int],
    detail_image_width: int,
    detail_image_height: int
    ) -> Dict[int, Dict[str, Any]]:
    """Fetches details for display. Now requires image dimensions."""
    if not cids: return {}
    logger.info(f"Fetching details for top {len(cids)} compound(s)...")
    compound_data = {cid: {} for cid in cids} # Pre-initialize

    # Properties Fetch
    properties_to_fetch = ["Title", "IUPACName", "IsomericSMILES", "InChI", "InChIKey"]
    cids_str = ','.join(map(str, cids))
    props_url = f"{BASE_URL}/compound/cid/{cids_str}/property/{','.join(properties_to_fetch)}/JSON"
    props_data = run_pubchem_search(props_url)
    if props_data and 'PropertyTable' in props_data and 'Properties' in props_data['PropertyTable']:
        for prop_entry in props_data['PropertyTable']['Properties']:
            cid = prop_entry.get('CID')
            if cid in compound_data: compound_data[cid].update(prop_entry)
    elif '_error_status' in props_data:
        logger.warning(f"Could not fetch properties: {props_data.get('_error_message')}")
    else:
        logger.warning("Could not parse properties response.")

    # Description Fetch
    desc_url = f"{BASE_URL}/compound/cid/{cids_str}/description/JSON"
    desc_data = run_pubchem_search(desc_url)
    if desc_data and 'InformationList' in desc_data and 'Information' in desc_data['InformationList']:
        valid_descriptions = [info for info in desc_data['InformationList']['Information'] if info.get('Description')]
        if valid_descriptions:
            for info in valid_descriptions:
                cid = info.get('CID')
                if cid in compound_data: compound_data[cid]['Description'] = info.get('Description')
    elif '_error_status' in desc_data:
        logger.warning(f"Could not fetch descriptions: {desc_data.get('_error_message')}")
    # else: logger.warning("Could not parse descriptions response.") # Less verbose


    # Add Links & Predict NMR for display
    logger.info("Adding links and predicting NMR for detailed display...")
    for cid in cids:
        if cid in compound_data:
            compound_data[cid]['PubChemLink'] = f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}"
            compound_data[cid]['ImageLink'] = f"https://pubchem.ncbi.nlm.nih.gov/image/imagefly.cgi?cid={cid}&width={detail_image_width}&height={detail_image_height}"
            # Re-fetch SMILES here if not guaranteed by property fetch? No, rely on props_data.
            smiles = compound_data[cid].get('IsomericSMILES')
            if RDKIT_AVAILABLE and smiles:
                 signals, integrations = get_hydrogen_environments(smiles)
                 compound_data[cid]['PredictedNMR'] = format_nmr_prediction(signals, integrations)
            else:
                 compound_data[cid]['PredictedNMR'] = "N/A" if RDKIT_AVAILABLE else "N/A (RDKit Missing)"

    logger.info("Detail fetching complete.")
    return compound_data