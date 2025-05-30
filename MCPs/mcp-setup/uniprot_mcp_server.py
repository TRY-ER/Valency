import json
import requests
from requests.adapters import HTTPAdapter, Retry
from fastmcp import FastMCP
import os
import sys
from env_loader import load_env_vars

load_env_vars()

# Create an MCP server
uniprot_host = os.getenv("UNIPROT_HOST", "0.0.0.0")
uniprot_port = int(os.getenv("UNIPROT_PORT", "8060")) # Assuming 8052 as a default, adjust if needed

mcp = FastMCP(
    "UniProt MCP Server",
    description="MCP server for querying the UniProt database."
)

UNIPROT_API_BASE_URL = "https://rest.uniprot.org"

# Configure requests session with retries
session = requests.Session()
retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[429, 500, 502, 503, 504])
session.mount("https://", HTTPAdapter(max_retries=retries))

@mcp.tool()
def search_uniprotkb(
    query_string: str,
    result_format: str = "json",
    fields: str = "",  # Changed default from None to ""
    size: int = 30,
    cursor: str = "",  # Changed default from None to ""
    include_isoform: bool = False,
    # compressed: bool = False # The requests library handles decompression automatically
) -> str:
    """
    Search the UniProtKB database. Look very cautiously at the query_string parameter for give valid queries.
    Args:
        query_string: The UniProt query string.
                      Examples:
                      - Simple term search: 'insulin'
                      - Multiple terms (AND): 'human antigen' or 'human AND antigen'
                      - Exact phrase: '"serum albumin"'
                      - Exclusion (NOT): 'human NOT mouse' or 'human -mouse'
                      - OR search: 'human OR mouse'
                      - Grouping: '(human OR mouse) AND kinase'
                      - Wildcard: 'anti*' (terms starting with 'anti'), '*kinase*' (terms containing 'kinase')
                      - Field-specific searches:
                          - Accession: 'accession:P62988'
                          - Active status: 'active:false' (for obsolete entries)
                          - Literature author: 'lit_author:ashburner'
                          - Protein name: 'protein_name:CD233'
                          - ChEBI ID: 'chebi:18420'
                          - PDB cross-reference count: 'xref_count_pdb:[20 TO *]'
                          - Creation date: 'date_created:[2012-10-01 TO *]'
                          - Modification date (active entries): 'date_modified:[2019-01-01 TO 2019-03-01] AND active:true'
                          - Sequence modification date: 'date_sequence_modified:[2012-01-01 TO 2012-03-01]'
                          - Database cross-reference: 'database:pfam' or 'xref:pdb-1aut'
                          - EC number: 'ec:3.2.1.23'
                          - Protein existence: 'existence:3'
                          - Protein family: 'family:serpin'
                          - Fragment status: 'fragment:true'
                          - Gene name: 'gene:HPSE' (matches HPSE, HPSE2)
                          - Exact gene name: 'gene_exact:HPSE' (matches only HPSE)
                          - Gene Ontology (GO) term: 'go:0015629' (actin cytoskeleton)
                          - Virus host ID: 'virus_host_id:10090' (mouse)
                          - InChIKey: 'inchikey:WQZGKKKJIJFFOK-GASJEMHNSA-N' (D-glucopyranose)
                          - Interactor: 'interactor:P00520'
                          - Keyword: 'keyword:toxin' or 'keyword:KW-0800'
                          - Sequence length: 'length:[500 TO 700]' or 'length:[50 TO *]'
                          - Sequence mass: 'mass:[500000 TO *]'
                          - Mass spectrometry method: 'cc_mass_spectrometry:maldi'
                          - Organelle: 'organelle:Mitochondrion'
                          - Organism name/ID: 'organism_name:"Ovis aries"' or 'organism_id:9940'
                          - Plasmid: 'plasmid:ColE1'
                          - Proteome ID: 'proteome:UP000005640' (human proteome)
                          - Proteome component: 'proteomecomponent:"chromosome 1" AND organism_id:9606'
                          - Secondary accession: 'sec_acc:P02023'
                          - Reviewed status (Swiss-Prot): 'reviewed:true'
                          - Publication scope: 'scope:mutagenesis'
                          - Strain: 'strain:wistar'
                          - Taxonomy name/ID: 'taxonomy_name:mammal' or 'taxonomy_id:40674'
                          - Tissue: 'tissue:liver'
                          - Web resource: 'cc_webresource:wikipedia'
                      - Escaping special characters: 'gene:L\\\\(1\\\\)2CB' (searches for gene name 'L(1)2CB')
        result_format: Desired format ('json', 'tsv', 'fasta', 'xml', 'txt', 'list', 'gff', 'obo', 'rdf', 'xlsx'). Defaults to 'json'.
        fields: Comma-separated list of column names to retrieve (applies to tsv, xlsx, json). E.g., 'id,xref_pdb,gene_names'.
        size: Number of results to retrieve per page (max 30 recommended). Defaults to 500.
        cursor: Cursor for pagination to retrieve the next page of results.
        include_isoform: Whether to include isoforms in the search results. Defaults to False.
    Returns:
        A JSON string containing the API response (data or error).
    """
    params = {"query": query_string, "format": result_format, "size": size}
    if fields:
        params["fields"] = fields
    if cursor:
        params["cursor"] = cursor
    if include_isoform:
        params["includeIsoform"] = "true" # API expects string "true"

    search_url = f"{UNIPROT_API_BASE_URL}/uniprotkb/search"
    response = None # Initialize response
    try:
        response = session.get(search_url, params=params)
        response.raise_for_status() # Raises an HTTPError for bad responses (4XX or 5XX)
        
        # For JSON, parse and re-dump to ensure consistent JSON string output
        if result_format == "json":
            return json.dumps(response.json())
        # For other text-based formats, return the text directly
        return json.dumps({"data": response.text, "content_type": response.headers.get("Content-Type")})

    except requests.exceptions.HTTPError as http_err:
        return json.dumps({"error": "HTTP error occurred", "details": str(http_err), "status_code": http_err.response.status_code, "response_text": http_err.response.text})
    except requests.exceptions.RequestException as req_err:
        return json.dumps({"error": "Request failed", "details": str(req_err)})
    except json.JSONDecodeError as json_err:
        # This might happen if result_format is json but the response isn't valid JSON (should be rare with UniProt)
        raw_response_text = response.text if response is not None else "Response object is None"
        return json.dumps({"error": "Failed to decode JSON response from UniProt", "details": str(json_err), "raw_response": raw_response_text})
    except Exception as e:
        return json.dumps({"error": "An unexpected error occurred", "details": str(e)})

@mcp.tool()
def get_uniprotkb_entry(uniprot_id: str, result_format: str = "json") -> str:
    """
    Retrieve a specific UniProtKB entry by its ID.
    Args:
        uniprot_id: The UniProtKB ID (e.g., 'P12345', 'SPIKE_SARS2').
        result_format: Desired format ('json', 'fasta', 'txt', 'xml', 'rdf', 'gff'). Defaults to 'json'.
    Returns:
        A JSON string containing the API response (data or error).
    """
    entry_url = f"{UNIPROT_API_BASE_URL}/uniprotkb/{uniprot_id}"
    params = {"format": result_format}
    response = None # Initialize response
    try:
        response = session.get(entry_url, params=params)
        response.raise_for_status()

        if result_format == "json":
            return json.dumps(response.json())
        return json.dumps({"data": response.text, "content_type": response.headers.get("Content-Type")})

    except requests.exceptions.HTTPError as http_err:
        return json.dumps({"error": "HTTP error occurred", "details": str(http_err), "status_code": http_err.response.status_code, "response_text": http_err.response.text})
    except requests.exceptions.RequestException as req_err:
        return json.dumps({"error": "Request failed", "details": str(req_err)})
    except json.JSONDecodeError as json_err:
        raw_response_text = response.text if response is not None else "Response object is None"
        return json.dumps({"error": "Failed to decode JSON response from UniProt", "details": str(json_err), "raw_response": raw_response_text})
    except Exception as e:
        return json.dumps({"error": "An unexpected error occurred", "details": str(e)})

if __name__ == "__main__":
    print(f"Starting UniProt MCP Server...")
    print(f"Server Name: {mcp.name}")
    sys.stdout.flush()
    mcp.run(transport="sse", host=uniprot_host, port=uniprot_port)
