import os
import json
import requests
import sys
from env_loader import load_env_vars

load_env_vars()

from fastmcp import FastMCP # Changed import

# --- Environment Variables for ADMET MCP Server ---
# Host and port for this MCP server
admet_mcp_host = os.getenv("ADMET_MCP_HOST", "0.0.0.0")
admet_mcp_port = int(os.getenv("ADMET_MCP_PORT", "8055")) # Port for uvicorn

# Host and port for the ADMET FastAPI application
admet_api_server_host = os.getenv("ADMET_SERVER_HOST", "0.0.0.0")
admet_api_server_port = int(os.getenv("ADMET_SERVER_PORT", "7373"))
admet_api_base_url = f"http://{admet_api_server_host}:{admet_api_server_port}"

mcp = FastMCP(
    "ADMET MCP Server",
    description="MCP server for retrieving ADMET predictions for a given SMILES string.",
    settings={"host": admet_mcp_host, "port": admet_mcp_port} # Configure host and port
    # The 'dependencies' argument might not be directly supported in jlowin/fastmcp's constructor in the same way.
    # If needed, refer to jlowin/fastmcp documentation for managing tool dependencies.
)

# --- ADMET API Tool ---

@mcp.tool()
def get_admet_prediction(smiles: str) -> str:
    """Get ADMET (Absorption, Distribution, Metabolism, Excretion, Toxicity)
    predictions for a given SMILES string.
    Args:
        smiles: The SMILES string representing the molecule.
    Returns:
        A JSON string containing the ADMET predictions, or an error message.
    """
    if not smiles or not isinstance(smiles, str):
        return json.dumps({"error": "Invalid input. SMILES string must be a non-empty string."})

    predict_url = f"{admet_api_base_url}/predict/"
    try:
        payload = {"smiles": smiles}
        headers = {"Content-Type": "application/json"} # Ensure correct content type
        
        # Using a timeout for the request
        response = requests.post(predict_url, json=payload, headers=headers, timeout=60) # 60 seconds timeout
        response.raise_for_status()  # Raises an HTTPError for bad responses (4XX or 5XX)
        
        return json.dumps(response.json())
        
    except requests.exceptions.HTTPError as http_err:
        error_detail = "No response content"
        if http_err.response is not None:
            try:
                error_detail = http_err.response.json() # Try to get JSON error from FastAPI
            except json.JSONDecodeError:
                error_detail = http_err.response.text # Fallback to text
        return json.dumps({"error": f"HTTP error occurred: {http_err}", "details": error_detail})
    except requests.exceptions.ConnectionError as conn_err:
        return json.dumps({"error": f"Connection error: Could not connect to ADMET service at {predict_url}.", "details": str(conn_err)})
    except requests.exceptions.Timeout as timeout_err:
        return json.dumps({"error": f"Request timed out while contacting ADMET service at {predict_url}.", "details": str(timeout_err)})
    except Exception as e:
        return json.dumps({"error": f"Failed to get ADMET prediction for SMILES '{smiles}'.", "details": str(e)})

# --- Main execution for direct run or mcp dev ---
if __name__ == "__main__":
    print(f"Starting ADMET MCP Server...")
    print(f"Server Name: {mcp.name}")
    print(f"ADMET FastAPI app expected at: {admet_api_base_url}")
    # The host and port for SSE are now passed to mcp.run()
    print(f"Attempting to run ADMET MCP Server with FastMCP SSE transport on host {admet_mcp_host}, port {admet_mcp_port}")
    sys.stdout.flush()

    # mcp.run() from jlowin/fastmcp will start a FastAPI server.
    # For SSE, specify transport, host, and port directly.
    mcp.run(transport="sse", host=admet_mcp_host, port=admet_mcp_port)

