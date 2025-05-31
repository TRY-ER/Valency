import os
import json
import requests
import sys 
from env_loader import load_env_vars

load_env_vars()

from fastmcp import FastMCP # Changed import

# Create an MCP server
alphafold_host = os.getenv("ALPHAFOLD_HOST", "0.0.0.0")
alphafold_port = int(os.getenv("ALPHAFOLD_PORT", "8056")) 
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
def get_alphafold_prediction(qualifier: str) -> str:
    """Get all AlphaFold models for a UniProt accession.
    Args:
        qualifier: UniProt accession (e.g., 'Q5VSL9').
    Returns:
        A JSON string containing the AlphaFold models data, or an error message.

    Example Succesful Returns Schema:
     [
        {
            "entryId": "string",
            "gene": "string",
            "sequenceChecksum": "string",
            "sequenceVersionDate": "string",
            "uniprotAccession": "string",
            "uniprotId": "string",
            "uniprotDescription": "string",
            "taxId": 0,
            "organismScientificName": "string",
            "uniprotStart": 0,
            "uniprotEnd": 0,
            "uniprotSequence": "string",
            "modelCreatedDate": "1941-09-22",
            "latestVersion": 0,
            "allVersions": [
            0
            ],
            "bcifUrl": "string",
            "cifUrl": "string",
            "pdbUrl": "string",
            "paeImageUrl": "string",
            "paeDocUrl": "string",
            "amAnnotationsUrl": "string",
            "amAnnotationsHg19Url": "string",
            "amAnnotationsHg38Url": "string",
            "isReviewed": true,
            "isReferenceProteome": true
        }
    ]  

    Example Value Error Returns Schema:

    {
        "detail": [
            {
            "loc": [
                "string",
                0
            ],
            "msg": "string",
            "type": "string"
            }
        ]
    }
    """
    response = None  # Initialize response
    try:
        url = f"{alphafold_api_base_url}/prediction/{qualifier}"
        params = {}
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

    Example Successful Returns Schema:
    {
        "uniprot_entry": {
            "ac": "string",
            "id": "string",
            "uniprot_checksum": "string",
            "sequence_length": 0,
            "segment_start": 0,
            "segment_end": 0
        },
        "structures": [
            {
            "summary": {
                "model_identifier": "string",
                "model_category": "EXPERIMENTALLY DETERMINED",
                "model_url": "string",
                "model_format": "PDB",
                "model_type": "ATOMIC",
                "model_page_url": "string",
                "provider": "string",
                "number_of_conformers": 0,
                "ensemble_sample_url": "string",
                "ensemble_sample_format": "PDB",
                "created": "string",
                "sequence_identity": 0,
                "uniprot_start": 0,
                "uniprot_end": 0,
                "coverage": 0,
                "experimental_method": "ELECTRON CRYSTALLOGRAPHY",
                "resolution": 0,
                "confidence_type": "pLDDT",
                "confidence_version": "string",
                "confidence_avg_local_score": 0,
                "oligomeric_state": "MONOMER",
                "preferred_assembly_id": "string",
                "entities": [
                {
                    "entity_type": "BRANCHED",
                    "entity_poly_type": "CYCLIC-PSEUDO-PEPTIDE",
                    "identifier": "string",
                    "identifier_category": "UNIPROT",
                    "description": "string",
                    "chain_ids": [
                    "string"
                    ]
                }
                ]
            }
            }
        ]
    }

    Example Value Error Returns Schema:
    {
        "detail": [
            {
            "loc": [
                "string",
                0
            ],
            "msg": "string",
            "type": "string"
            }
        ]
    }
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

    Example Successful Returns Schema:
    {
        "accession": "string",
        "id": "string",
        "sequence": "string",
        "annotation": [
            {
            "type": "MUTAGEN",
            "description": "string",
            "source_name": "string",
            "source_url": "string",
            "evidence": "COMPUTATIONAL/PREDICTED",
            "residues": [
                0
            ],
            "regions": [
                {
                "start": 1,
                "end": 120,
                "annotation_value": [
                    "string"
                ],
                "unit": "string"
                }
            ]
            }
        ]
    }

    Example Value Error Returns Schema:
    {
        "detail": [
            {
            "loc": [
                "string",
                0
            ],
            "msg": "string",
            "type": "string"
            }
       ]
    }
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
