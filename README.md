# PubChem Isomer Finder & Analyzer

This web application allows users to search for chemical compounds (isomers) based on molecular formula, specific structural features (functional groups), and predicted 1H NMR characteristics (number of signals and integration ratios). It leverages the PubChem PUG REST API for compound data and (optionally) RDKit for cheminformatics calculations.

## Core Chemical Concepts

### 1. Isomerism
Isomers are molecules that have the same molecular formula but different structural or spatial arrangements of atoms. This application primarily focuses on **structural isomers**, where atoms are connected in different orders.

*   **Molecular Formula:** The fundamental starting point. For example, C₄H₁₀ can represent *n*-butane or isobutane.
*   **Structural Constraints (Functional Groups):** Users can specify required functional groups (e.g., "must contain a carboxylic acid," "must have an aromatic ring"). This significantly narrows down the search space from all possible isomers of a given formula to only those possessing the desired chemical moieties. The constraints are defined using SMARTS (SMiles ARbitrary Target Specification), a language for describing molecular patterns.
*   **1H NMR Spectroscopy Data:**
    *   **Number of Signals:** In a 1H NMR spectrum, chemically non-equivalent protons give rise to distinct signals. The number of signals provides information about the symmetry of the molecule.
    *   **Integration Ratio:** The area under each signal (integration) is proportional to the number of protons contributing to that signal. The simplified ratio of these integrations (e.g., 3:2:1 for 9 protons split into groups of 3, 2, and 1 equivalent protons, simplified from 6H:4H:2H) is a powerful constraint for distinguishing isomers.

### 2. SMARTS (SMiles ARbitrary Target Specification)
SMARTS is used to define the structural constraints. It's a powerful query language that allows for the specification of substructures, atom types, bond types, and environmental properties. For example:
*   `[CX3](=O)[OX2H1]` represents a carboxylic acid group.
*   `c1ccccc1` represents a benzene ring.
The application uses a predefined library of common functional groups mapped to their SMARTS patterns.

### 3. RDKit for NMR Prediction
If RDKit is available, the application uses it to:
1.  **Parse SMILES:** Convert SMILES strings (obtained from PubChem) into RDKit molecule objects.
2.  **Add Hydrogens:** Explicitly add hydrogen atoms to the molecular graph, as they are often implicit in SMILES.
3.  **Determine Proton Equivalency:** RDKit's `CanonicalRankAtoms` function is used to assign a rank to each atom in the molecule. Atoms (including hydrogens) with the same rank are considered chemically equivalent in terms of NMR.
4.  **Count Signals:** The number of unique ranks among hydrogen atoms corresponds to the predicted number of 1H NMR signals.
5.  **Calculate Integration:** For each unique hydrogen rank, the number of hydrogens sharing that rank is counted. This provides the raw integration values.
6.  **Simplify Ratio:** The raw integration values are simplified by dividing by their greatest common divisor (GCD) to produce the characteristic integration ratio.

This prediction is a simplified model and doesn't account for complex phenomena like diastereotopic protons splitting or magnetic anisotropy effects that can further complicate real spectra. However, it provides a good first-order approximation for structural elucidation.

## Application Architecture & Design

The application is a Flask-based web server with a distinct front-end (HTML/CSS/JavaScript) and back-end (Python/Flask).

### Backend (Python/Flask - `app.py` and `core_logic` package)

1.  **Flask App (`app.py`):**
    *   **Routing:** Defines URL endpoints for the main search page (`/`), results display (`/results`), Server-Sent Events stream (`/search-stream`), and AJAX detail fetching (`/fetch-batch-details`).
    *   **Request Handling:** Parses user input from forms and AJAX requests.
    *   **Context Processing:** Injects global variables (like RDKit availability) into templates.
    *   **Server-Sent Events (SSE) for Search Progress:**
        *   The `/search-stream` route handles the multi-step search process.
        *   It uses a generator function with `@stream_with_context` to send progress updates to the client in real-time without keeping a single HTTP request open indefinitely or resorting to polling.
        *   Messages (log, status, warning, error, completion) are formatted as JSON and sent over the event stream.
        *   A `core_logic_progress_callback` is passed to long-running backend functions to queue messages for the SSE stream.
    *   **Interaction with `core_logic`:** Delegates chemical search and data processing tasks to the `core_logic` package.
    *   **Templating:** Uses Jinja2 to render HTML pages, passing data retrieved from `core_logic`.
    *   **Custom Jinja Filter:** A `fromjson` filter is registered to allow parsing JSON strings directly within templates.

2.  **Core Logic (`core_logic` package):**
    This package encapsulates all the chemical data handling and API interactions. It's structured into several modules for better organization:

    *   **`__init__.py`:**
        *   Handles RDKit import and sets `RDKIT_AVAILABLE` flag.
        *   Loads the `SMARTS_LIBRARY_DATA` from `data/smarts_library.json`.
        *   Provides a `get_smarts_categories` helper function.
        *   Exports key functions and constants from submodules for easier access.

    *   **`data/smarts_library.json`:**
        *   Stores the extensive library of functional groups, their names, SMARTS patterns, categories, and descriptions in a JSON format, making it easier to maintain and extend than embedding it directly in Python code.

    *   **`pubchem_api.py`:**
        *   `run_pubchem_search()`: A robust wrapper around `requests` to interact with the PubChem PUG REST API.
        *   Handles GET and POST requests, JSON parsing, error checking (HTTP errors, JSON decode errors), and implements a retry mechanism with exponential backoff for transient server-side issues.
        *   Returns a dictionary, embedding error information (`_error_status`, `_error_message`) if an error occurs.
        *   Defines constants like `BASE_URL`, `REQUEST_TIMEOUT`.

    *   **`search_operations.py`:**
        *   `find_candidates_by_structure()`:
            1.  Performs an initial PubChem search by molecular formula to get a `CacheKey` (a temporary ID for a list of results on PubChem servers).
            2.  If structural constraints (SMARTS) are provided, it iteratively refines this list on PubChem's servers by applying each SMARTS pattern as a substructure search against the current `CacheKey`. For each step except the last, it requests a new `CacheKey`; for the last step, it requests the final list of CIDs.
            3.  Communicates progress via the `progress_callback`.

    *   **`nmr_processing.py`:**
        *   `get_hydrogen_environments()`: (RDKit dependent)
            *   Takes a SMILES string.
            *   Uses RDKit to parse SMILES, add explicit hydrogens.
            *   Calls `Chem.CanonicalRankAtoms(mol_h, breakTies=False, includeChirality=True, includeIsotopes=False)` to determine atom equivalency. The `breakTies=False` is crucial for differentiating protons that might be equivalent under simpler ranking schemes but are distinct in NMR (e.g., geminal protons in a chiral environment if chirality isn't considered, or homotopic vs enantiotopic vs diastereotopic protons). `includeChirality=True` helps in distinguishing diastereotopic protons. `includeIsotopes=False` is used as we are generally interested in the common 1H.
            *   Counts unique ranks among hydrogen atoms to get the number of signals.
            *   Counts hydrogens per unique rank for integration values.
        *   `calculate_predicted_ratio()`: Calculates the simplified proton ratio from integration counts by finding the GCD of the integration values and dividing each by it.
        *   `format_nmr_prediction()`: Formats the NMR data into a human-readable string.
        *   `filter_candidates_by_nmr()`: (RDKit dependent)
            1.  Fetches SMILES strings in batches for the candidate CIDs (using `fetch_smiles_batch`).
            2.  For each compound with a valid SMILES:
                *   Predicts NMR signals and ratio using `get_hydrogen_environments` and `calculate_predicted_ratio`.
                *   Filters the compound list based on user-provided `required_signals` and `required_ratio`.
            3.  Reports progress and warnings.

    *   **`compound_utils.py`:**
        *   `fetch_smiles_batch()`:
            *   Efficiently fetches SMILES strings (Isomeric preferred, fallback to Canonical) for a list of CIDs from PubChem in batches to avoid overly long URLs or request timeouts.
            *   Handles API errors and missing SMILES data.
        *   `fetch_compound_details()`:
            *   Fetches detailed information (Title, IUPAC Name, SMILES, InChIKey, Molecular Formula, Molecular Weight, Description, image link) for a given list of CIDs from PubChem.
            *   If RDKit is available, it also calls `get_hydrogen_environments` and `format_nmr_prediction` to add predicted NMR data to the details.
            *   Imports NMR functions locally within the function to resolve circular dependencies with `nmr_processing.py`.

### Frontend (HTML/CSS/JavaScript - `templates` and `static` directories)

1.  **Templates (`templates/`):**
    *   `layout.html`: Base template providing common HTML structure, Bootstrap CSS/JS, and custom stylesheet. Defines blocks for title, content, and extra scripts/styles.
    *   `index.html`: The main search page.
        *   Renders the form with inputs for formula, structural constraints (dynamically generated from `SMARTS_CATEGORIES`), and NMR constraints.
        *   Includes JavaScript to handle form submission via SSE.
        *   Displays a progress bar and a log area for real-time updates from the server during the search.
        *   Upon search completion (via SSE 'completion' event), it shows a "Show Results" button instead of auto-redirecting, using the URL provided in the SSE payload.
    *   `results.html`: Displays the search results.
        *   Shows a summary of the search parameters used.
        *   Displays any warnings or errors.
        *   Initially renders a few detailed result cards (using the `_result_card.html` partial).
        *   Implements a "Load More" button if more results are available than initially shown. Clicking this button triggers an AJAX call to `/fetch-batch-details` to get data for the next set of CIDs.
        *   The JavaScript dynamically creates new result cards from the AJAX response and appends them to the page.
        *   Displays a summary list of all found CIDs, with each CID hyperlinked to its PubChem page. This list is truncated for display if very long.
    *   `_result_card.html`: A Jinja2 partial template responsible for rendering a single compound's detail card. This promotes reusability between the initial page load and AJAX-loaded results.

2.  **Static Files (`static/`):**
    *   `css/style.css`: Custom CSS rules for styling beyond Bootstrap defaults.

3.  **Client-Side JavaScript (embedded in `index.html` and `results.html`):**
    *   **`index.html` JS:**
        *   Handles form submission, preventing default browser action.
        *   Initiates an `EventSource` connection to `/search-stream`.
        *   Listens for `onopen`, `onmessage`, and custom `completion` events from the SSE stream.
        *   Updates the progress bar and log list based on received messages.
        *   On `completion` (type `results_ready`): enables a "Show Results" button and stores the redirect URL. If (type `error_final`): displays the error and allows retrying.
        *   Handles `onerror` for the EventSource.
    *   **`results.html` JS:**
        *   Handles the "Load More" functionality.
        *   Takes the full list of CIDs (passed as JSON from the server).
        *   When "Load More" is clicked, it slices the next batch of CIDs from the full list.
        *   Makes an AJAX `POST` request to `/fetch-batch-details` with the batch of CIDs.
        *   On success, receives JSON data for the new compounds.
        *   Dynamically calls `createResultCardHtml()` (a JS function) to build HTML for each new result card and appends it to the DOM.
        *   `createResultCardHtml()`: A utility function to construct the HTML string for a result card, ensuring consistent rendering. It includes basic HTML escaping for security.

## Key Design Decisions & Trade-offs

*   **SSE for Progress:** Chosen for its efficiency in providing real-time, unidirectional updates from server to client over a single HTTP connection, suitable for long-running search tasks.
*   **`core_logic` Package:** Separates chemical logic from web framework concerns, improving modularity and testability.
*   **JSON for SMARTS Library:** Makes the SMARTS data easier to manage and extend than hardcoding a large Python dictionary.
*   **Client-Side "Load More":** The server sends all found CIDs to the client initially. The client then requests details for batches on demand.
    *   *Pro:* Simpler initial backend logic for the results page; reduces initial data transfer for details.
    *   *Con:* If an extremely large number of CIDs (e.g., >10,000s) are found, sending all CIDs as a JSON string might be slow or hit URL/data limits if it were passed in GET. (Currently, it's passed to the template context, which is fine). The main limitation would be browser memory if the JS array of CIDs becomes excessively large.
*   **RDKit as Optional Dependency:** The application functions for formula and structural search even without RDKit, but NMR features are disabled. This increases portability.
*   **Error Handling:** Attempts to provide informative error messages both in server logs and to the client via SSE or the results page.
*   **Delayed Import for Circular Dependency:** `compound_utils.py` imports `nmr_processing` functions locally within `fetch_compound_details` to resolve a circular import, a common pattern in Python.
