import os
import requests
import json
from mcp.server.fastmcp import FastMCP
from typing import Optional
from dotenv import load_dotenv

load_dotenv()

DEEP_PK_API_URL = "https://biosig.lab.uq.edu.au/deeppk/api/predict"
DEFAULT_MCP_HOST = "0.0.0.0"
DEFAULT_MCP_PORT = 8053

# Define constants for MCP server metadata
MCP_TITLE = "Deep-PK ADMET Prediction MCP Server"
MCP_DESCRIPTION = "Provides tools to predict ADMET properties of drug molecules using the Deep-PK API."
MCP_DEPENDENCIES = ["requests"]

mcp_host = os.getenv("DEEP_PK_MCP_HOST", DEFAULT_MCP_HOST)
mcp_port = int(os.getenv("DEEP_PK_MCP_PORT", str(DEFAULT_MCP_PORT)))

mcp = FastMCP(
    MCP_TITLE,
    description=MCP_DESCRIPTION,
    dependencies=MCP_DEPENDENCIES,
    host=mcp_host,
    port=mcp_port
)

@mcp.tool()
def submit_deep_pk_job(
    input_molecule: str,
    input_type: str,
    email: Optional[str] = None,
    pred_type: str = "ADMET"
) -> str:
    """Submit a prediction job to the Deep-PK API.
    Provide input molecule via the 'input_molecule' argument and specify its type using 'input_type'.
    'input_type' can be one of: 'smiles_string', 'smiles_file_path', or 'sdf_file_path'.
    File paths must be accessible to the server running this MCP.

    Args:
        input_molecule (str): The molecule input. This can be a SMILES string,
                              a path to a file containing SMILES strings, or a path to an SDF file.
        input_type (str): Specifies the type of input_molecule.
                          Options: "smiles_string", "smiles_file_path", "sdf_file_path".
        email (str, optional): Email for contact when the job is finished.
        pred_type (str, optional): Type of prediction. Options: "Absorption", "Distribution",
                                   "Metabolism", "Excretion", "ADMET" (default).

    Returns:
        str: A JSON string containing the job_id if successful, or an error message.
    """
    multipart_payload = {}  # All form fields will go here to ensure multipart/form-data
    opened_files = []  # To keep track of files to close

    if input_type == "smiles_string":
        multipart_payload['smiles'] = (None, input_molecule)
    elif input_type == "smiles_file_path":
        try:
            file_handle = open(input_molecule, 'rb')
            opened_files.append(file_handle)
            multipart_payload['smiles_file'] = (os.path.basename(input_molecule), file_handle)
        except IOError as e:
            return json.dumps({"error": f"Failed to open smiles_file at '{input_molecule}': {e}"})
    elif input_type == "sdf_file_path":
        try:
            file_handle = open(input_molecule, 'rb')
            opened_files.append(file_handle)
            multipart_payload['sdf_file'] = (os.path.basename(input_molecule), file_handle)
        except IOError as e:
            return json.dumps({"error": f"Failed to open sdf_file at '{input_molecule}': {e}"})
    else:
        for f in opened_files:
            f.close()
        return json.dumps({"error": "Invalid input_type. Choose from: 'smiles_string', 'smiles_file_path', 'sdf_file_path'."})

    if email:
        multipart_payload['email'] = (None, email)

    valid_pred_types = ["Absorption", "Distribution", "Metabolism", "Excretion", "ADMET"]
    pred_type_for_api = pred_type.lower()
    valid_pred_types_lower = [v.lower() for v in valid_pred_types]

    if pred_type_for_api not in valid_pred_types_lower:
        for f in opened_files:
            f.close()
        return json.dumps({"error": f"Invalid pred_type '{pred_type}'. Choose from: {', '.join(valid_pred_types)}"})
    
    multipart_payload['pred_type'] = (None, pred_type_for_api)

    response_obj = None
    try:
        response_obj = requests.post(DEEP_PK_API_URL, files=multipart_payload)
        response_obj.raise_for_status()
        return json.dumps(response_obj.json())
    except requests.exceptions.HTTPError as e:
        err_text = e.response.text if e.response else "No response text"
        return json.dumps({"error": f"HTTP error occurred: {e} - {err_text}"})
    except json.JSONDecodeError as e:
        raw_text = response_obj.text if response_obj is not None else "No response object or text available"
        return json.dumps({"error": "Failed to decode JSON response from API.", "raw_text": raw_text, "details": str(e)})
    except requests.exceptions.RequestException as e:
        return json.dumps({"error": f"API request failed: {e}"})
    except IOError as e:
        return json.dumps({"error": f"I/O error during request: {e}"})
    except Exception as e:
        return json.dumps({"error": "An unexpected error occurred.", "details": str(e)})
    finally:
        for f in opened_files:
            f.close()

@mcp.tool()
def get_deep_pk_results(job_id: str) -> str:
    """Retrieve the status or results of a Deep-PK prediction job.

    Args:
        job_id (str): The unique ID of the job.

    Returns:
        str: A JSON string containing the job status or the prediction results if completed, or an error message.
    """
    params = {'job_id': job_id}
    response_obj = None
    try:
        response_obj = requests.get(DEEP_PK_API_URL, params=params)
        response_obj.raise_for_status()
        return json.dumps(response_obj.json())
    except requests.exceptions.HTTPError as e:
        err_text = e.response.text if e.response else "No response text"
        status_code = e.response.status_code if e.response else None
        return json.dumps({"error": f"HTTP error occurred: {e} - {err_text}", "status_code": status_code})
    except requests.exceptions.JSONDecodeError:
        raw_text = response_obj.text if response_obj is not None else "No response object or text available"
        return json.dumps({"error": "Failed to decode JSON response from API. The response was not valid JSON.", 
                "raw_text": raw_text})
    except requests.exceptions.RequestException as e:
        return json.dumps({"error": f"API request failed: {e}"})
    except Exception as e:
        return json.dumps({"error": "An unexpected error occurred while fetching results.", "details": str(e)})

if __name__ == "__main__":
    mcp_transport_type = os.getenv("MCP_TRANSPORT", "stdio")
    
    print(f"Starting Deep-PK MCP Server...")
    print(f"MCP Name: {MCP_TITLE}")
    print(f"Description: {MCP_DESCRIPTION}")
    print(f"Dependencies: {', '.join(MCP_DEPENDENCIES)}")
    
    if mcp_transport_type == "stdio":
        print(f"Transport: stdio")
        mcp.run(transport="stdio")
    elif mcp_transport_type == "sse":
        print(f"Transport: sse")
        print(f"Host: {mcp_host}")
        print(f"Port: {mcp_port}")
        print(f"Deep-PK MCP Server will be available via SSE.")
        mcp.run(transport="sse")
    else:
        print(f"Error: Unknown MCP_TRANSPORT type '{mcp_transport_type}'. Supported types are 'stdio' and 'sse'.")
        print("Defaulting to stdio transport.")
        mcp.run(transport="stdio")

