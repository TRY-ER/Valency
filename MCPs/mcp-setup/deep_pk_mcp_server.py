import os
import requests
import json
from mcp.server.fastmcp import FastMCP
from typing import Optional  # Added for Optional type hint
from dotenv import load_dotenv  # Import load_dotenv

load_dotenv()  # Load .env file at the beginning

# Configuration for the Deep-PK API and MCP Server
DEEP_PK_API_URL = "https://biosig.lab.uq.edu.au/deeppk/api/predict"
DEFAULT_MCP_HOST = "0.0.0.0"
DEFAULT_MCP_PORT = 8053  # Using a different port than the RCSB example

mcp_host = os.getenv("DEEP_PK_MCP_HOST", DEFAULT_MCP_HOST)
mcp_port = int(os.getenv("DEEP_PK_MCP_PORT", str(DEFAULT_MCP_PORT)))

mcp = FastMCP(
    "Deep-PK ADMET Prediction MCP Server",
    description="Provides tools to predict ADMET properties of drug molecules using the Deep-PK API.",
    dependencies=["requests"],
    host=mcp_host,
    port=mcp_port
)

@mcp.tool()
def submit_deep_pk_job(
    smiles_string: Optional[str] = None, 
    smiles_file_path: Optional[str] = None, 
    sdf_file_path: Optional[str] = None, 
    email: Optional[str] = None, 
    pred_type: str = "ADMET"
) -> str:  # Changed return type to str
    """Submit a prediction job to the Deep-PK API.
    Provide input molecule(s) via one of the following: smiles_string, smiles_file_path, or sdf_file_path.
    File paths must be accessible to the server running this MCP.

    Args:
        smiles_string (str, optional): Single SMILES string.
        smiles_file_path (str, optional): Path to a file containing SMILES strings.
        sdf_file_path (str, optional): Path to an SDF file.
        email (str, optional): Email for contact when the job is finished.
        pred_type (str, optional): Type of prediction. Options: "Absorption", "Distribution", 
                                   "Metabolism", "Excretion", "ADMET" (default).

    Returns:
        str: A JSON string containing the job_id if successful, or an error message.
    """
    files_to_send = {}
    data_payload = {}
    opened_files = []  # To keep track of files to close

    if smiles_string:
        data_payload['smiles'] = smiles_string
    elif smiles_file_path:
        try:
            # API expects 'smiles_file' as the field name for the file
            file_handle = open(smiles_file_path, 'rb')
            opened_files.append(file_handle)
            files_to_send['smiles_file'] = (os.path.basename(smiles_file_path), file_handle)
        except IOError as e:
            return json.dumps({"error": f"Failed to open smiles_file at '{smiles_file_path}': {e}"})
    elif sdf_file_path:
        try:
            # API expects 'sdf_file' as the field name for the file
            file_handle = open(sdf_file_path, 'rb')
            opened_files.append(file_handle)
            files_to_send['sdf_file'] = (os.path.basename(sdf_file_path), file_handle)
        except IOError as e:
            return json.dumps({"error": f"Failed to open sdf_file at '{sdf_file_path}': {e}"})
    else:
        return json.dumps({"error": "Input required: Provide one of smiles_string, smiles_file_path, or sdf_file_path."})

    if email:
        data_payload['email'] = email

    valid_pred_types = ["Absorption", "Distribution", "Metabolism", "Excretion", "ADMET"]
    if pred_type not in valid_pred_types:
        # Close any opened files before returning error
        for f in opened_files:
            f.close()
        return json.dumps({"error": f"Invalid pred_type '{pred_type}'. Choose from: {', '.join(valid_pred_types)}"})
    data_payload['pred_type'] = pred_type

    response_obj = None  # Initialize response_obj to None
    try:
        response_obj = requests.post(DEEP_PK_API_URL, data=data_payload, files=files_to_send if files_to_send else None)
        response_obj.raise_for_status()  # Raise an HTTPError for bad responses (4XX or 5XX)
        return json.dumps(response_obj.json())  # Return JSON string
    except IOError as e:  # Specific error for file opening
        return json.dumps({"error": f"File I/O error: {e}"})  # Return JSON string
    except requests.exceptions.HTTPError as e:
        err_text = e.response.text if e.response else "No response text"
        return json.dumps({"error": f"HTTP error occurred: {e} - {err_text}"})  # Return JSON string
    except requests.exceptions.JSONDecodeError:  # Catch JSONDecodeError specifically
        raw_text = response_obj.text if response_obj is not None else "No response object or text available"
        return json.dumps({"error": "Failed to decode JSON response from API.", "raw_text": raw_text})  # Return JSON string
    except requests.exceptions.RequestException as e:  # Catch other request-related errors
        return json.dumps({"error": f"API request failed: {e}"})  # Return JSON string
    except Exception as e:  # Generic catch-all for other unexpected errors
        return json.dumps({"error": "An unexpected error occurred.", "details": str(e)})  # Return JSON string
    finally:
        for f in opened_files:
            f.close()

@mcp.tool()
def get_deep_pk_results(job_id: str) -> str:  # Changed return type to str
    """Retrieve the status or results of a Deep-PK prediction job.

    Args:
        job_id (str): The unique ID of the job.

    Returns:
        str: A JSON string containing the job status or the prediction results if completed, or an error message.
    """
    params = {'job_id': job_id}
    response_obj = None  # Initialize response_obj to None
    try:
        response_obj = requests.get(DEEP_PK_API_URL, params=params)
        response_obj.raise_for_status()  # Raise an HTTPError for bad responses
        
        # The response content can be a JSON object for status or results.
        # Example status: {"status": "running"}
        # Example result: {"0": {"SMILES": "...", ...}}
        return json.dumps(response_obj.json())  # Return JSON string
    except requests.exceptions.HTTPError as e:
        err_text = e.response.text if e.response else "No response text"
        status_code = e.response.status_code if e.response else None
        return json.dumps({"error": f"HTTP error occurred: {e} - {err_text}", "status_code": status_code})  # Return JSON string
    except requests.exceptions.JSONDecodeError:  # Catch JSONDecodeError specifically and before RequestException
        raw_text = response_obj.text if response_obj is not None else "No response object or text available"
        return json.dumps({"error": "Failed to decode JSON response from API. The response was not valid JSON.", 
                "raw_text": raw_text})  # Return JSON string
    except requests.exceptions.RequestException as e:  # Catch other request-related errors
        return json.dumps({"error": f"API request failed: {e}"})  # Return JSON string
    except Exception as e:  # Generic catch-all for other unexpected errors
        return json.dumps({"error": "An unexpected error occurred while fetching results.", "details": str(e)})  # Return JSON string

if __name__ == "__main__":
    mcp_transport_type = os.getenv("MCP_TRANSPORT", "stdio")
    
    print(f"Starting Deep-PK MCP Server...")
    print(f"MCP Name: {mcp.name}")
    print(f"Description: {mcp.description}")
    print(f"Dependencies: {', '.join(mcp.dependencies)}")
    
    if mcp_transport_type == "stdio":
        print(f"Transport: stdio")
        mcp.run(transport="stdio")
    elif mcp_transport_type == "sse":
        print(f"Transport: sse")
        print(f"Host: {mcp.host}")
        print(f"Port: {mcp.port}")
        print(f"Deep-PK MCP Server will be available via SSE.")
        mcp.run(transport="sse")  # host and port are taken from mcp instance
    else:
        print(f"Error: Unknown MCP_TRANSPORT type '{mcp_transport_type}'. Supported types are 'stdio' and 'sse'.")
        print("Defaulting to stdio transport.")
        mcp.run(transport="stdio")

