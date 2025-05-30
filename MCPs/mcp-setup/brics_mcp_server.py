import os
import json
import sys
from env_loader import load_env_vars

load_env_vars()

from fastmcp import FastMCP
from utils.generators.BRICSGenerator import BRICSGenerator

# --- Environment Variables for BRICS MCP Server ---
brics_mcp_host = os.getenv("BRICS_MCP_HOST", "0.0.0.0")
brics_mcp_port = int(os.getenv("BRICS_MCP_PORT", "8058")) # New port for this server

mcp = FastMCP(
    "BRICS MCP Server",
    description="MCP server for generating molecular candidates using BRICS algorithm.",
    settings={"host": brics_mcp_host, "port": brics_mcp_port}
)

# --- BRICS Generation Tool ---

@mcp.tool()
def get_brics_candidates(smiles_list: list[str], is_polymer: bool = False) -> str:
    """Generate molecular candidates from a list of SMILES strings using the BRICS algorithm.
    Args:
        smiles_list: A list of SMILES strings.
        is_polymer: A boolean flag indicating if the input SMILES are polymers (default: False).
    Returns:
        A JSON string containing the generated candidates and their count, or an error message.

    Example Response Schema:
    {
        "candidates": [
            "C1=CC=CC=C1",
            "C1=CC=CC=C2C=CC=CC=C2"
        ],
        "count": 2
    }
    """
    if not smiles_list or not isinstance(smiles_list, list) or not all(isinstance(s, str) for s in smiles_list):
        return json.dumps({"error": "Invalid input. smiles_list must be a non-empty list of strings."})

    try:
        generator = BRICSGenerator(verbose=False) # verbose can be set based on server config if needed
        candidates, count = generator.generate(smiles_list, is_polymer=is_polymer)
        return json.dumps({"candidates": candidates, "count": count})
    except ValueError as ve:
        return json.dumps({"error": f"Input validation error: {str(ve)}"})
    except Exception as e:
        return json.dumps({"error": f"Failed to generate BRICS candidates.", "details": str(e)})

# --- Main execution for direct run or mcp dev ---
if __name__ == "__main__":
    print(f"Starting BRICS MCP Server...")
    print(f"Server Name: {mcp.name}")
    print(f"Attempting to run BRICS MCP Server with FastMCP SSE transport on host {brics_mcp_host}, port {brics_mcp_port}")
    sys.stdout.flush()

    mcp.run(transport="sse", host=brics_mcp_host, port=brics_mcp_port)
