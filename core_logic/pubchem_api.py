import requests
import json
import time
import logging
from typing import Optional, Dict, Any

logger = logging.getLogger(__name__)

# --- Constants ---
BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
REQUEST_TIMEOUT = 120  # seconds
SMILES_BATCH_SIZE = 500 # Max CIDs per SMILES request

# --- API Call Helper ---
def run_pubchem_search(
    url: str,
    params: Optional[Dict] = None,
    request_type: str = "GET",
    data: Optional[Dict] = None, # For POST requests
    retries: int = 2,
    delay: float = 2.0
) -> Dict[str, Any]: # Return type changed to always be a dict, error info embedded
    """
    Helper function to make requests to the PubChem PUG REST API.
    Returns a dictionary with the JSON response or an error structure.
    Error structure: {'_error_status': int, '_error_message': str}
    """
    headers = {'Accept': 'application/json', 'User-Agent': 'IsomerWebApp/1.3'} # Updated version
    last_exception = None
    last_status_code = None

    logger.debug(f"PubChem Request ({request_type}): {url} | Params: {params} | Data: {data}")

    for attempt in range(retries + 1):
        try:
            if request_type.upper() == "GET":
                response = requests.get(url, params=params, headers=headers, timeout=REQUEST_TIMEOUT)
            elif request_type.upper() == "POST":
                headers['Content-Type'] = 'application/x-www-form-urlencoded' # Common for PUG REST POST
                response = requests.post(url, data=data, headers=headers, timeout=REQUEST_TIMEOUT)
            else:
                logger.error(f"Unsupported request type '{request_type}' for URL: {url}")
                return {'_error_status': 400, '_error_message': f"Internal Error: Unsupported request type '{request_type}'"}

            last_status_code = response.status_code
            response.raise_for_status()  # Raises HTTPError for 4xx/5xx responses
            
            # Handle cases where response might be empty but successful (e.g., 204 No Content)
            if not response.content:
                 logger.debug(f"PubChem request successful with empty content (Status: {last_status_code}) for URL: {url}")
                 return {} # Or some other indicator of success with no data, if appropriate for context

            result_data = response.json()
            logger.debug(f"PubChem request successful (Status: {last_status_code}) for URL: {url}")
            return result_data

        except requests.exceptions.Timeout as e:
            last_exception = e
            logger.warning(f"PubChem request timed out (Attempt {attempt+1}/{retries+1}). URL: {url}")
        except requests.exceptions.HTTPError as e:
            last_exception = e
            status_code = response.status_code if response is not None else -1
            server_error_message = f"PubChem API Error {status_code}"
            
            try:
                error_details_json = response.json()
                fault_msg = error_details_json.get('Fault', {}).get('Message')
                details_list = error_details_json.get('Fault', {}).get('Details')
                
                details_str = None
                if isinstance(details_list, list):
                    details_str = '; '.join(details_list)
                elif isinstance(details_list, str):
                    details_str = details_list

                if fault_msg: server_error_message += f": {fault_msg}"
                if details_str: server_error_message += f" ({details_str})"
            except json.JSONDecodeError:
                logger.warning(f"Could not decode JSON error response from PubChem. Status: {status_code}. Response text: {response.text[:200] if response else 'N/A'}")
            except Exception as inner_ex:
                logger.warning(f"Error parsing PubChem error details (Status {status_code}): {inner_ex}")

            if 500 <= status_code <= 599 and attempt < retries:
                logger.warning(f"PubChem API returned HTTP {status_code} (Attempt {attempt+1}/{retries+1}). Retrying... URL: {url}")
            else:
                logger.error(f"PubChem API returned non-retryable HTTP Error {status_code} for URL {url} (Final Attempt). Message: {server_error_message}")
                return {'_error_status': status_code, '_error_message': server_error_message}
        
        except json.JSONDecodeError as e:
            last_exception = e
            logger.error(f"Failed to decode JSON response from PubChem for {url} (Attempt {attempt+1}/{retries+1}). Status: {last_status_code}. Response: {response.text[:200] if response else 'N/A'}")
            # Return error structure instead of None
            return {'_error_status': 500, '_error_message': f"Invalid JSON response from PubChem (expected success but got non-JSON). Response text: {response.text[:200] if response else 'N/A'}"}

        except requests.exceptions.RequestException as e:
            last_exception = e
            logger.error(f"Network or request error occurred: {e} (Attempt {attempt+1}/{retries+1}). URL: {url}")
        except Exception as e: # Catch any other unexpected errors
            last_exception = e
            logger.exception(f"An unexpected error occurred during API call to {url} (Attempt {attempt+1}/{retries+1}).")


        if attempt < retries:
            calculated_delay = delay * (1.5**attempt) # Exponential backoff
            logger.info(f"Waiting {calculated_delay:.2f}s before retry for {url}...")
            time.sleep(calculated_delay)
        # Implicitly, if loop finishes without returning, it means all retries failed

    # Fallback after all retries failed
    error_message = f"Request failed after {retries+1} attempts. Last status: {last_status_code}. Last error: {type(last_exception).__name__ if last_exception else 'N/A'}"
    if last_exception and hasattr(last_exception, 'args') and last_exception.args:
        error_message += f" ({str(last_exception.args[0])[:100]})"
    logger.error(f"{error_message} URL: {url}")
    return {'_error_status': last_status_code or 500, '_error_message': error_message}