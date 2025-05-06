import json
import os
import logging

# Setup Logging for the core_logic package
logger = logging.getLogger(__name__)

# --- RDKit Import ---
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem # Though not explicitly used in provided snippet, good to have if needed
    RDKIT_AVAILABLE = True
    logger.info("RDKit successfully imported.")
except ImportError:
    RDKIT_AVAILABLE = False
    logger.warning("RDKit not found. NMR-related features will be disabled.")

# --- Load SMARTS Library Data ---
_smarts_library_path = os.path.join(os.path.dirname(__file__), 'data', 'smarts_library.json')
try:
    with open(_smarts_library_path, 'r') as f:
        SMARTS_LIBRARY_DATA = json.load(f)
    logger.info(f"Successfully loaded SMARTS library from {_smarts_library_path}")
except FileNotFoundError:
    logger.critical(f"SMARTS library file not found at: {_smarts_library_path}")
    SMARTS_LIBRARY_DATA = {}
except json.JSONDecodeError as e:
    logger.critical(f"Error decoding SMARTS library JSON from {_smarts_library_path}: {e}")
    SMARTS_LIBRARY_DATA = {}
except Exception as e:
    logger.critical(f"An unexpected error occurred while loading SMARTS library: {e}")
    SMARTS_LIBRARY_DATA = {}

# --- Helper to get categories (can be used by app.py) ---
def get_smarts_categories(library_data=None):
    if library_data is None:
        library_data = SMARTS_LIBRARY_DATA

    categories = {}
    for key, data in library_data.items():
        cat = data.get('category', 'Uncategorized')
        if cat not in categories:
            categories[cat] = []
        categories[cat].append({
            'key': key,
            'name': data.get('name', key),
            'desc': data.get('desc', '')
        })
    # Sort items within each category
    for cat in categories:
        categories[cat].sort(key=lambda x: x['name'])
    # Sort categories themselves
    return dict(sorted(categories.items()))


# --- Expose necessary functions and variables ---
# For cleaner imports: from core_logic import find_candidates_by_structure, etc.
from .pubchem_api import run_pubchem_search, BASE_URL, REQUEST_TIMEOUT, SMILES_BATCH_SIZE
from .search_operations import find_candidates_by_structure
from .nmr_processing import (
    get_hydrogen_environments,
    calculate_predicted_ratio,
    format_nmr_prediction,
    filter_candidates_by_nmr
)
from .compound_utils import fetch_smiles_batch, fetch_compound_details

__all__ = [
    'RDKIT_AVAILABLE', 'SMARTS_LIBRARY_DATA', 'get_smarts_categories',
    'run_pubchem_search', 'BASE_URL', 'REQUEST_TIMEOUT', 'SMILES_BATCH_SIZE',
    'find_candidates_by_structure',
    'get_hydrogen_environments', 'calculate_predicted_ratio', 'format_nmr_prediction',
    'filter_candidates_by_nmr',
    'fetch_smiles_batch', 'fetch_compound_details'
]