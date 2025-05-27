import os
import json
import requests
import sys 
from env_loader import load_env_vars

load_env_vars()

from fastmcp import FastMCP # Changed import

# Create an MCP server
alphafold_host = os.getenv("ALPHAFOLD_HOST", "0.0.0.0")
alphafold_port = int(os.getenv("ALPHAFOLD_PORT", "8052")) 
alphafold_api_base_url = os.getenv("ALPHAFOLD_API_BASE_URL", "https://alphafold.ebi.ac.uk/api")

mcp = FastMCP(
    "AlphaFold MCP Server", # Name of the server
    description="MCP server for interacting with the AlphaFold Protein Structure Database API.",
    settings={"host": alphafold_host, "port": alphafold_port} # Configure host and port
    # The 'dependencies' argument might not be directly supported in jlowin/fastmcp's constructor in the same way.
    # If needed, refer to jlowin/fastmcp documentation for managing tool dependencies.
)

# --- AlphaFold API Tools ---

@mcp.tool()
def get_alphafold_prediction(qualifier: str, sequence_checksum: str | None = None) -> str:
    """Get all AlphaFold models for a UniProt accession.
    Args:
        qualifier: UniProt accession (e.g., 'Q5VSL9').
        sequence_checksum: Optional CRC64 checksum of the UniProt sequence.
    Returns:
        A JSON string containing the AlphaFold models data, or an error message.
    """
    response = None  # Initialize response
    try:
        url = f"{alphafold_api_base_url}/prediction/{qualifier}"
        params = {}
        if sequence_checksum:
            params["sequence_checksum"] = sequence_checksum
        
        response = requests.get(url, params=params)
        response.raise_for_status()  # Raises an HTTPError for bad responses (4XX or 5XX)
        return json.dumps(response.json())
    except requests.exceptions.HTTPError as http_err:
        return json.dumps({"error": f"HTTP error occurred: {http_err}", "details": response.text if response else "No response"})
    except Exception as e:
        return json.dumps({"error": f"Failed to get AlphaFold prediction for qualifier '{qualifier}'.", "details": str(e)})

@mcp.tool()
def get_uniprot_summary(qualifier: str) -> str:
    """Get summary details for a UniProt residue range.
    Args:
        qualifier: UniProtKB accession number (AC), entry name (ID), or CRC64 checksum of the UniProt sequence (e.g., 'Q5VSL9').
    Returns:
        A JSON string containing summary details for the UniProt residue range, or an error message.
    """
    response = None  # Initialize response
    try:
        url = f"{alphafold_api_base_url}/uniprot/summary/{qualifier}.json"
        response = requests.get(url)
        response.raise_for_status()
        return json.dumps(response.json())
    except requests.exceptions.HTTPError as http_err:
        return json.dumps({"error": f"HTTP error occurred: {http_err}", "details": response.text if response else "No response"})
    except Exception as e:
        return json.dumps({"error": f"Failed to get UniProt summary for qualifier '{qualifier}'.", "details": str(e)})

@mcp.tool()
def get_alphafold_annotations(qualifier: str, annotation_type: str) -> str:
    """Get all annotations for a UniProt residue range.
    Args:
        qualifier: UniProt accession (e.g., 'Q5VSL9').
        annotation_type: Type of annotation (e.g., 'MUTAGEN' for AlphaMissense).
    Returns:
        A JSON string containing annotations for the UniProt residue range, or an error message.
    """
    response = None  # Initialize response
    try:
        url = f"{alphafold_api_base_url}/annotations/{qualifier}"
        params = {"annotation_type": annotation_type}
        response = requests.get(url, params=params)
        response.raise_for_status()
        return json.dumps(response.json())
    except requests.exceptions.HTTPError as http_err:
        return json.dumps({"error": f"HTTP error occurred: {http_err}", "details": response.text if response else "No response"})
    except Exception as e:
        return json.dumps({"error": f"Failed to get AlphaFold annotations for qualifier '{qualifier}' with type '{annotation_type}'.", "details": str(e)})

# --- Main execution for direct run or mcp dev ---
if __name__ == "__main__":
    print(f"Starting AlphaFold MCP Server...")
    print(f"Server Name: {mcp.name}")
    # The host and port for SSE are now passed to mcp.run()
    print(f"Attempting to run AlphaFold MCP Server with FastMCP SSE transport on host {alphafold_host}, port {alphafold_port}")
    sys.stdout.flush()

    # mcp.run() from jlowin/fastmcp will start a FastAPI server.
    # For SSE, specify transport, host, and port directly.
    mcp.run(transport="sse", host=alphafold_host, port=alphafold_port)
