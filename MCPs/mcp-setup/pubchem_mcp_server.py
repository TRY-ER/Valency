import os
import json
import requests # Dependency for making HTTP requests
import urllib.parse # For URL encoding
import sys 

from env_loader import load_env_vars 
from fastmcp import FastMCP # Changed import

load_env_vars() 

# --- Server Setup ---
PUBCHEM_BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

pubchem_mcp_host = os.getenv("PUBCHEM_MCP_HOST", "0.0.0.0")
pubchem_mcp_port = int(os.getenv("PUBCHEM_MCP_PORT", "8054")) 

mcp = FastMCP(
    "PubChem MCP Server",
    # Removed dependencies from here, ensure it's installed in the environment if needed by tools
    # Removed host/port from settings, will be passed to run()
    description="MCP server for querying the PubChem PUG REST API."
)

# --- Helper Function for API Calls ---
def _fetch_pubchem_data(endpoint: str, params: dict | None = None) -> str:
    """Helper function to fetch data from PubChem PUG REST API.
    Returns a JSON string with a standardized structure:
    Success: {"status": "success", "result": <data_payload>}
    Error:   {"status": "error", "error": {"code": "...", "message": "...", "details": ...}}
    """
    url = f"{PUBCHEM_BASE_URL}{endpoint}"
    try:
        response = requests.get(url, params=params, timeout=30)
        response.raise_for_status()  # Raises an HTTPError for bad responses (4XX or 5XX)
        
        content_type = response.headers.get("Content-Type", "").lower()

        if "application/json" in content_type:
            return json.dumps({
                "status": "success",
                "result": response.json()
            })
        elif "text/plain" in content_type:
             return json.dumps({
                "status": "success",
                "result": {"text_content": response.text, "original_content_type": response.headers.get("Content-Type", "")}
            })
        else: # Handle other content types (e.g., PNG, XML if not directly requested as such by a specialized tool)
            return json.dumps({
                "status": "error",
                "error": {
                    "code": "UNEXPECTED_CONTENT_TYPE",
                    "message": f'PubChem API returned an unexpected content type: {response.headers.get("Content-Type", "N/A")}.',
                    "details": {"response_text_preview": response.text[:200] if response.text else None}
                }
            })

    except requests.exceptions.HTTPError as http_err:
        error_details_payload = None
        status_code_str = "Unknown"
        if http_err.response is not None:
            status_code_str = str(http_err.response.status_code)
            try:
                error_details_payload = http_err.response.json()
            except json.JSONDecodeError:
                error_details_payload = http_err.response.text
        return json.dumps({
            "status": "error",
            "error": {
                "code": "PUBCHEM_API_ERROR",
                "message": f"PubChem API HTTP error: {status_code_str}",
                "details": error_details_payload
            }
        })
    except requests.exceptions.Timeout:
        return json.dumps({
            "status": "error",
            "error": {
                "code": "TIMEOUT_ERROR",
                "message": "PubChem API request timed out."
            }
        })
    except requests.exceptions.RequestException as req_err:
        return json.dumps({
            "status": "error",
            "error": {
                "code": "REQUEST_ERROR",
                "message": "PubChem API request failed.",
                "details": str(req_err)
            }
        })
    except Exception as e: # Catch any other unexpected errors within this helper
        return json.dumps({
            "status": "error",
            "error": {
                "code": "MCP_HELPER_ERROR",
                "message": "An unexpected error occurred within the PubChem MCP helper function.",
                "details": str(e)
            }
        })

# --- Compound Tools ---

@mcp.tool()
def get_compound_by_cid(cid: str) -> str:
    """Retrieve the full compound record from PubChem by its Compound ID (CID).
    Args:
        cid: The PubChem Compound ID (e.g., '2244' for aspirin).
    Returns:
        A JSON string with a status and either the compound data in 'result' or error details in 'error'.
    """
    if not cid or not isinstance(cid, str) or not cid.strip():
        return json.dumps({
            "status": "error",
            "error": {"code": "VALIDATION_ERROR", "message": "CID must be a non-empty string."}
        })
    endpoint = f"/compound/cid/{urllib.parse.quote(cid.strip())}/record/JSON"
    return _fetch_pubchem_data(endpoint)

@mcp.tool()
def get_cids_by_name(name: str, name_type: str | None = None) -> str:
    """Search for PubChem Compound IDs (CIDs) by chemical name.
    Args:
        name: The chemical name to search for (e.g., 'aspirin').
        name_type: Optional. Set to 'word' to match words within the name. Default is exact match.
    Returns:
        A JSON string with status and a list of CIDs in 'result' or error details.
    """
    if not name or not isinstance(name, str) or not name.strip():
        return json.dumps({
            "status": "error",
            "error": {"code": "VALIDATION_ERROR", "message": "Name must be a non-empty string."}
        })
    if name_type and not isinstance(name_type, str):
        return json.dumps({
            "status": "error",
            "error": {"code": "VALIDATION_ERROR", "message": "name_type must be a string if provided."}
        })
        
    safe_name = urllib.parse.quote(name.strip())
    endpoint = f"/compound/name/{safe_name}/cids/JSON"
    params = {}
    if name_type:
        params["name_type"] = name_type.strip()
    return _fetch_pubchem_data(endpoint, params=params if params else None)

@mcp.tool()
def get_compound_properties(cids: list[str], properties_list: list[str]) -> str:
    """Retrieve specified properties for one or more PubChem CIDs.
    Args:
        cids: A list of PubChem Compound IDs (e.g., ['2244', '31200']).
        properties_list: A list of properties to retrieve (e.g., ['MolecularWeight', 'InChIKey']).
    Returns:
        A JSON string with status and properties in 'result' or error details.
    """
    if not cids or not isinstance(cids, list) or not all(isinstance(cid, str) and cid.strip() for cid in cids):
        return json.dumps({
            "status": "error",
            "error": {"code": "VALIDATION_ERROR", "message": "cids must be a non-empty list of non-empty strings."}
        })
    if not properties_list or not isinstance(properties_list, list) or not all(isinstance(p, str) and p.strip() for p in properties_list):
        return json.dumps({
            "status": "error",
            "error": {"code": "VALIDATION_ERROR", "message": "properties_list must be a non-empty list of non-empty strings."}
        })

    safe_cids = ",".join([urllib.parse.quote(cid.strip()) for cid in cids])
    safe_properties = ",".join([urllib.parse.quote(prop.strip()) for prop in properties_list])
    
    endpoint = f"/compound/cid/{safe_cids}/property/{safe_properties}/JSON"
    return _fetch_pubchem_data(endpoint)

@mcp.tool()
def get_compound_synonyms_by_cid(cid: str) -> str:
    """Retrieve all synonyms for a given PubChem Compound ID (CID).
    Args:
        cid: The PubChem Compound ID (e.g., '2244').
    Returns:
        A JSON string with status and synonyms in 'result' or error details.
    """
    if not cid or not isinstance(cid, str) or not cid.strip():
        return json.dumps({
            "status": "error",
            "error": {"code": "VALIDATION_ERROR", "message": "CID must be a non-empty string."}
        })
    endpoint = f"/compound/cid/{urllib.parse.quote(cid.strip())}/synonyms/JSON"
    return _fetch_pubchem_data(endpoint)

@mcp.tool()
def get_compound_image_pubchem_url(input_type: str, identifier: str) -> str:
    """Constructs and returns the direct PubChem URL for a compound's PNG image.
    Args:
        input_type: The type of identifier. Examples: 'cid', 'name', 'smiles', 'inchikey'.
        identifier: The actual identifier string.
    Returns:
        A JSON string with status and the image URL in 'result.image_url' or error details.
    """
    if not input_type or not isinstance(input_type, str) or not input_type.strip():
        return json.dumps({
            "status": "error",
            "error": {"code": "VALIDATION_ERROR", "message": "Input type must be a non-empty string."}
        })
    if not identifier or not isinstance(identifier, str) or not identifier.strip():
        return json.dumps({
            "status": "error",
            "error": {"code": "VALIDATION_ERROR", "message": "Identifier must be a non-empty string."}
        })

    input_type_clean = input_type.strip().lower()
    identifier_clean = identifier.strip()
    valid_input_types = ['cid', 'name', 'smiles', 'inchikey']

    if input_type_clean not in valid_input_types:
        return json.dumps({
            "status": "error",
            "error": {
                "code": "VALIDATION_ERROR",
                "message": f"Invalid input_type '{input_type_clean}'. Must be one of {valid_input_types}."
            }
        })

    safe_identifier = urllib.parse.quote_plus(identifier_clean)
    image_url = f"{PUBCHEM_BASE_URL}/compound/{input_type_clean}/{safe_identifier}/PNG"
    
    return json.dumps({
        "status": "success",
        "result": {"image_url": image_url}
    })

@mcp.tool()
def get_cids_by_smiles(smiles: str) -> str:
    """Retrieve PubChem CIDs for a given SMILES string.
    Args:
        smiles: The SMILES string (e.g., 'CCO' for ethanol).
    Returns:
        A JSON string with status and a list of CIDs in 'result' or error details.
    """
    if not smiles or not isinstance(smiles, str) or not smiles.strip():
        return json.dumps({
            "status": "error",
            "error": {"code": "VALIDATION_ERROR", "message": "SMILES string must be non-empty."}
        })
    safe_smiles = urllib.parse.quote(smiles.strip())
    endpoint = f"/compound/smiles/{safe_smiles}/cids/JSON"
    return _fetch_pubchem_data(endpoint)

@mcp.tool()
def get_cids_by_inchikey(inchikey: str) -> str:
    """Retrieve PubChem CIDs for a given InChIKey.
    Args:
        inchikey: The InChIKey (e.g., 'LFQSCWFLJHTTHZ-UHFFFAOYSA-N' for ethanol).
    Returns:
        A JSON string with status and a list of CIDs in 'result' or error details.
    """
    if not inchikey or not isinstance(inchikey, str) or not inchikey.strip():
        return json.dumps({
            "status": "error",
            "error": {"code": "VALIDATION_ERROR", "message": "InChIKey must be a non-empty string."}
        })
    safe_inchikey = urllib.parse.quote(inchikey.strip())
    endpoint = f"/compound/inchikey/{safe_inchikey}/cids/JSON"
    return _fetch_pubchem_data(endpoint)

@mcp.tool()
def fast_identity_search_by_cid(cid: str, identity_type: str = "same_connectivity") -> str:
    """Perform a fast identity search for compounds related to a given CID.
    Args:
        cid: The PubChem Compound ID (e.g., '2244').
        identity_type: Type of identity search. Default is 'same_connectivity'.
    Returns:
        A JSON string with status and a list of CIDs in 'result' or error details.
    """
    if not cid or not isinstance(cid, str) or not cid.strip():
        return json.dumps({
            "status": "error",
            "error": {"code": "VALIDATION_ERROR", "message": "CID must be a non-empty string."}
        })
    if not identity_type or not isinstance(identity_type, str) or not identity_type.strip():
        return json.dumps({
            "status": "error",
            "error": {"code": "VALIDATION_ERROR", "message": "identity_type must be a non-empty string."}
        })

    valid_identity_types = [
        "same_connectivity", "same_stereo", "same_isotope", "same_tautomer", 
        "same_stereo_isotope", "same_connectivity_stereo_isotope", "similar_stereo_isotope"
    ]
    identity_type_clean = identity_type.strip().lower()
    if identity_type_clean not in valid_identity_types:
        return json.dumps({
            "status": "error",
            "error": {"code": "VALIDATION_ERROR", "message": f"Invalid identity_type. Choose from {valid_identity_types}"}
        })
        
    endpoint = f"/compound/fastidentity/cid/{urllib.parse.quote(cid.strip())}/cids/JSON"
    params = {"identity_type": identity_type_clean}
    return _fetch_pubchem_data(endpoint, params=params)

@mcp.tool()
def fast_substructure_search_by_smiles(smiles: str, strip_hydrogen: bool | None = None) -> str:
    """Perform a fast substructure search using a SMILES string.
    Args:
        smiles: The SMILES string for the substructure query.
        strip_hydrogen: Optional. If True, hydrogens are stripped.
    Returns:
        A JSON string with status and matching CIDs in 'result' or error details.
    """
    if not smiles or not isinstance(smiles, str) or not smiles.strip():
        return json.dumps({
            "status": "error",
            "error": {"code": "VALIDATION_ERROR", "message": "SMILES string must be non-empty."}
        })
    if strip_hydrogen is not None and not isinstance(strip_hydrogen, bool):
         return json.dumps({
            "status": "error",
            "error": {"code": "VALIDATION_ERROR", "message": "strip_hydrogen must be a boolean if provided."}
        })

    safe_smiles = urllib.parse.quote(smiles.strip())
    endpoint = f"/compound/fastsubstructure/smiles/{safe_smiles}/cids/JSON"
    params = {}
    if strip_hydrogen is not None:
        params["StripHydrogen"] = str(strip_hydrogen).lower()
    return _fetch_pubchem_data(endpoint, params=params if params else None)

@mcp.tool()
def fast_similarity_2d_search_by_cid(cid: str, threshold: int = 90) -> str:
    """Perform a fast 2D similarity search based on a CID.
    Args:
        cid: The PubChem CID to find similar compounds to.
        threshold: Tanimoto similarity threshold (0-100). Default is 90.
    Returns:
        A JSON string with status and similar CIDs in 'result' or error details.
    """
    if not cid or not isinstance(cid, str) or not cid.strip():
        return json.dumps({
            "status": "error",
            "error": {"code": "VALIDATION_ERROR", "message": "CID must be a non-empty string."}
        })
    if not isinstance(threshold, int) or not (0 <= threshold <= 100):
        return json.dumps({
            "status": "error",
            "error": {"code": "VALIDATION_ERROR", "message": "Threshold must be an integer between 0 and 100."}
        })

    endpoint = f"/compound/fastsimilarity_2d/cid/{urllib.parse.quote(cid.strip())}/cids/JSON"
    params = {"Threshold": threshold}
    return _fetch_pubchem_data(endpoint, params=params)

@mcp.tool()
def get_cids_by_xref(xref_type: str, xref_value: str) -> str:
    """Retrieve PubChem CIDs by a cross-reference (XRef).
    Args:
        xref_type: The type of cross-reference (e.g., 'PatentID', 'RegistryID').
        xref_value: The value of the cross-reference.
    Returns:
        A JSON string with status and CIDs in 'result' or error details.
    """
    if not xref_type or not isinstance(xref_type, str) or not xref_type.strip():
        return json.dumps({
            "status": "error",
            "error": {"code": "VALIDATION_ERROR", "message": "XRef type must be a non-empty string."}
        })
    if not xref_value or not isinstance(xref_value, str) or not xref_value.strip():
        return json.dumps({
            "status": "error",
            "error": {"code": "VALIDATION_ERROR", "message": "XRef value must be a non-empty string."}
        })
    
    safe_xref_type = urllib.parse.quote(xref_type.strip())
    safe_xref_value = urllib.parse.quote(xref_value.strip())
    
    endpoint = f"/compound/xref/{safe_xref_type}/{safe_xref_value}/cids/JSON"
    return _fetch_pubchem_data(endpoint)

@mcp.tool()
def get_cids_by_mass(
    mass_type: str, 
    value_or_min: float, 
    max_value: float | None = None
    # mass_unit is not directly supported in URL path by PUG REST for these lookups
) -> str:
    """Retrieve PubChem CIDs by mass.
    Args:
        mass_type: Type of mass: 'molecular_weight', 'exact_mass', or 'monoisotopic_mass'.
        value_or_min: The exact mass value or the minimum mass for a range query.
        max_value: Optional. The maximum mass for a range query.
    Returns:
        A JSON string with status and CIDs in 'result' or error details.
    """
    if not mass_type or not isinstance(mass_type, str) or not mass_type.strip():
        return json.dumps({
            "status": "error",
            "error": {"code": "VALIDATION_ERROR", "message": "mass_type must be a non-empty string."}
        })
    if not isinstance(value_or_min, (float, int)):
        return json.dumps({
            "status": "error",
            "error": {"code": "VALIDATION_ERROR", "message": "value_or_min must be a number."}
        })
    if max_value is not None and not isinstance(max_value, (float, int)):
        return json.dumps({
            "status": "error",
            "error": {"code": "VALIDATION_ERROR", "message": "max_value must be a number if provided."}
        })

    valid_mass_types = ["molecular_weight", "exact_mass", "monoisotopic_mass"]
    mass_type_clean = mass_type.strip().lower()
    if mass_type_clean not in valid_mass_types:
        return json.dumps({
            "status": "error",
            "error": {"code": "VALIDATION_ERROR", "message": f"Invalid mass_type. Choose from {valid_mass_types}."}
        })

    safe_mass_type = urllib.parse.quote(mass_type_clean)
    
    if max_value is not None: # Range query
        if value_or_min > max_value:
            return json.dumps({
                "status": "error",
                "error": {"code": "VALIDATION_ERROR", "message": "value_or_min cannot be greater than max_value for a range query."}
            })
        endpoint = f"/compound/{safe_mass_type}/range/{urllib.parse.quote(str(value_or_min))}/{urllib.parse.quote(str(max_value))}/cids/JSON"
    else: # Exact match query (using range with min=max as per previous logic)
        endpoint = f"/compound/{safe_mass_type}/range/{urllib.parse.quote(str(value_or_min))}/{urllib.parse.quote(str(value_or_min))}/cids/JSON"
        
    return _fetch_pubchem_data(endpoint)

# --- Main execution for direct run or mcp dev ---
if __name__ == "__main__":
    # transport = os.getenv("MCP_TRANSPORT", "sse") # No longer needed for direct SSE run

    print(f"Starting PubChem MCP Server...")
    print(f"Server Name: {mcp.name}")
    sys.stdout.flush()

    # Default to SSE transport with host and port passed directly
    print(f"Attempting to run PubChem MCP Server with FastMCP SSE transport on host {pubchem_mcp_host}, port {pubchem_mcp_port}")
    sys.stdout.flush()
    mcp.run(transport="sse", host=pubchem_mcp_host, port=pubchem_mcp_port)

