import requests
from urllib.parse import quote

# Step 1: Search for compounds with the molecular formula C3H6O2 and retrieve a cachekey
formula_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastformula/C3H6O2/cids/JSON?list_return=cachekey"
response = requests.get(formula_url)

if response.status_code != 200:
    print(f"Error in formula search: {response.status_code}")
    print(response.text)
    exit()

try:
    data = response.json()
    cachekey = data['IdentifierList']['CacheKey']
    print(f"CacheKey obtained: {cachekey}")
except KeyError:
    print("CacheKey not found in response.")
    print(response.text)
    exit()

# Step 2: Substructure search for the ester functional group using SMARTS pattern
substructure_smarts = "[OX2]C(=O)[#6]"  # SMARTS pattern for ester group
encoded_smarts = quote(substructure_smarts)
substructure_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsubstructure/smarts/{encoded_smarts}/cids/JSON?cachekey={cachekey}"

response_sub = requests.get(substructure_url)

if response_sub.status_code == 200:
    data_sub = response_sub.json()
    cids = data_sub.get('IdentifierList', {}).get('CID', [])
    if cids:
        print(f"Found {len(cids)} ester compounds with formula C3H6O2:")
        for cid in cids:
            print(f"CID: {cid}")
    else:
        print("No esters found with the formula C3H6O2.")
else:
    print(f"Error in substructure search: {response_sub.status_code}")
    print(response_sub.text)