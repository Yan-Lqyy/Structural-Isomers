import logging
from typing import List, Dict, Optional, Tuple, Callable, Any

from .pubchem_api import run_pubchem_search, BASE_URL, SMILES_BATCH_SIZE
# from .nmr_processing import get_hydrogen_environments, format_nmr_prediction # MOVED
from . import RDKIT_AVAILABLE

logger = logging.getLogger(__name__)

# --- Batch SMILES Fetching ---
def fetch_smiles_batch(
    cids: List[int],
    progress_callback: Optional[Callable[[str], None]] = None
) -> Tuple[Dict[int, str], Optional[str]]: # Returns (cid_smiles_map, warning_message_or_none)
    """
    Fetches Isomeric (preferred) or Canonical SMILES for a list of CIDs in batches.
    """
    def report_progress(message: str):
        logger.info(f"[SMILESFetch] {message}")
        if progress_callback:
            try:
                progress_callback(message)
            except Exception as e:
                logger.warning(f"Progress callback failed for message '{message[:50]}...': {e}")

    if not cids:
        return {}, None

    report_progress(f"Initiating SMILES fetch for {len(cids)} CIDs (batch size: {SMILES_BATCH_SIZE})...")
    
    cid_smiles_map: Dict[int, str] = {}
    all_warnings: List[str] = []
    
    total_batches = (len(cids) + SMILES_BATCH_SIZE - 1) // SMILES_BATCH_SIZE

    for i in range(0, len(cids), SMILES_BATCH_SIZE):
        batch_num = (i // SMILES_BATCH_SIZE) + 1
        batch_cids = cids[i : i + SMILES_BATCH_SIZE]
        
        if batch_num == 1 or batch_num == total_batches or batch_num % (max(1, total_batches // 10)) == 0 :
            report_progress(f"  Fetching SMILES batch {batch_num}/{total_batches} ({len(batch_cids)} CIDs)...")

        cids_str = ','.join(map(str, batch_cids))
        properties_to_fetch = "IsomericSMILES,CanonicalSMILES"
        url = f"{BASE_URL}/compound/cid/{cids_str}/property/{properties_to_fetch}/JSON"
        
        smiles_data = run_pubchem_search(url)

        if '_error_status' in smiles_data:
            error_msg = smiles_data.get('_error_message', f'API error fetching SMILES for batch {batch_num}')
            logger.error(error_msg)
            all_warnings.append(f"Failed to fetch SMILES for CIDs in batch {batch_num}: {error_msg[:100]}...")
            continue 

        if not smiles_data or 'PropertyTable' not in smiles_data or 'Properties' not in smiles_data['PropertyTable']:
            warn_msg = f"Unexpected response structure for SMILES batch {batch_num}. Data: {str(smiles_data)[:150]}..."
            logger.warning(warn_msg)
            all_warnings.append(warn_msg)
            continue

        props_found_in_batch = 0
        for prop_entry in smiles_data['PropertyTable']['Properties']:
            cid_val = prop_entry.get('CID')
            if cid_val is None: 
                logger.warning("Found property entry without CID in SMILES batch.")
                continue
            
            if cid_val not in batch_cids:
                logger.warning(f"Received SMILES for unexpected CID {cid_val} in batch {batch_num}.")
                continue

            iso_smiles = prop_entry.get('IsomericSMILES')
            can_smiles = prop_entry.get('CanonicalSMILES')
            
            valid_smiles = None
            if iso_smiles and iso_smiles.strip() and iso_smiles.lower() != 'n/a':
                valid_smiles = iso_smiles
            elif can_smiles and can_smiles.strip() and can_smiles.lower() != 'n/a':
                valid_smiles = can_smiles
            
            if valid_smiles:
                cid_smiles_map[cid_val] = valid_smiles
                props_found_in_batch += 1
            else:
                logger.debug(f"No valid SMILES (Isomeric or Canonical) found for CID {cid_val} in batch {batch_num}.")
        
        if props_found_in_batch < len(batch_cids):
            missing_count = len(batch_cids) - props_found_in_batch
            all_warnings.append(f"Missing valid SMILES for {missing_count} CIDs in batch {batch_num}.")

    num_retrieved = len(cid_smiles_map)
    num_requested = len(cids)
    if num_retrieved < num_requested:
        all_warnings.append(f"Could not retrieve valid SMILES for {num_requested - num_retrieved} out of {num_requested} CIDs.")
    
    report_progress(f"SMILES fetching complete. Retrieved valid SMILES for {num_retrieved}/{num_requested} CIDs.")
    
    final_warning_message = "; ".join(all_warnings) if all_warnings else None
    return cid_smiles_map, final_warning_message


# --- Detail Fetching ---
def fetch_compound_details(
    cids: List[int],
    detail_image_width: int,
    detail_image_height: int
) -> Dict[int, Dict[str, Any]]:
    """
    Fetches detailed information for a list of CIDs.
    Includes properties, description, image link, and predicted NMR (if RDKit available).
    """
    # Import here to break circular dependency
    from .nmr_processing import get_hydrogen_environments, format_nmr_prediction

    if not cids:
        return {}

    logger.info(f"Fetching details for {len(cids)} compound(s)...")
    compound_data: Dict[int, Dict[str, Any]] = {cid: {'CID': cid} for cid in cids} 

    cids_str = ','.join(map(str, cids))

    # 1. Fetch Properties
    properties_to_fetch = ["Title", "IUPACName", "IsomericSMILES", "CanonicalSMILES", "InChI", "InChIKey", "MolecularFormula", "MolecularWeight"]
    props_url = f"{BASE_URL}/compound/cid/{cids_str}/property/{','.join(properties_to_fetch)}/JSON"
    props_response = run_pubchem_search(props_url)

    if '_error_status' in props_response:
        logger.warning(f"Could not fetch properties for CIDs: {props_response.get('_error_message')}")
    elif props_response and 'PropertyTable' in props_response and 'Properties' in props_response['PropertyTable']:
        props_found_count = 0
        for prop_entry in props_response['PropertyTable']['Properties']:
            cid_val = prop_entry.get('CID')
            if cid_val in compound_data:
                compound_data[cid_val].update(prop_entry)
                props_found_count +=1
        logger.info(f"Fetched properties for {props_found_count}/{len(cids)} CIDs.")
    else:
        logger.warning("Malformed response or no properties found for CIDs.")

    # 2. Fetch Descriptions
    desc_url = f"{BASE_URL}/compound/cid/{cids_str}/description/JSON"
    desc_response = run_pubchem_search(desc_url)
    
    if '_error_status' in desc_response:
        logger.warning(f"Could not fetch descriptions: {desc_response.get('_error_message')}")
    elif desc_response and 'InformationList' in desc_response and 'Information' in desc_response['InformationList']:
        desc_found_count = 0
        processed_cids_for_desc = set() 
        for info_entry in desc_response['InformationList']['Information']:
            cid_val = info_entry.get('CID')
            if cid_val in compound_data and cid_val not in processed_cids_for_desc:
                description = info_entry.get('Description')
                source_name = info_entry.get('DescriptionSourceName') 
                if description:
                    compound_data[cid_val]['Description'] = description
                    if source_name:
                         compound_data[cid_val]['DescriptionSourceName'] = source_name
                    processed_cids_for_desc.add(cid_val)
                    desc_found_count +=1
        logger.info(f"Fetched descriptions for {desc_found_count}/{len(cids)} CIDs.")
    else:
        logger.warning("Malformed response or no descriptions found.")
        
    # 3. Add Links and Predict NMR
    logger.info("Augmenting details with links and NMR predictions...")
    nmr_predictions_successful = 0
    for cid_val, details in compound_data.items(): 
        details['PubChemLink'] = f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid_val}"
        details['ImageLink'] = f"https://pubchem.ncbi.nlm.nih.gov/image/imagefly.cgi?cid={cid_val}&width={detail_image_width}&height={detail_image_height}"

        smiles_for_nmr = details.get('IsomericSMILES')
        if not smiles_for_nmr or smiles_for_nmr.lower() == 'n/a':
            smiles_for_nmr = details.get('CanonicalSMILES')

        if RDKIT_AVAILABLE and smiles_for_nmr and smiles_for_nmr.lower() != 'n/a':
            signals, integrations = get_hydrogen_environments(smiles_for_nmr)
            details['PredictedNMR'] = format_nmr_prediction(signals, integrations)
            if signals is not None: 
                nmr_predictions_successful +=1
        elif not RDKIT_AVAILABLE:
            details['PredictedNMR'] = "N/A (RDKit Missing)"
        else: 
            details['PredictedNMR'] = "N/A (SMILES unavailable for NMR)"
            
    logger.info(f"NMR predictions processed for {nmr_predictions_successful}/{len(cids)} compounds in details.")
    logger.info("Compound detail fetching and processing complete.")
    return compound_data