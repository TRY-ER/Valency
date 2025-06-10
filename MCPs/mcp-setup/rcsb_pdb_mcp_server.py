import os
from env_loader import load_env_vars # Changed from dotenv to env_loader

load_env_vars() # Load variables from .env file using custom loader

from fastmcp import FastMCP # Changed import
from rcsbsearchapi.search import (
    TextQuery,
    AttributeQuery,
    SequenceQuery,
    SeqMotifQuery,
    StructSimilarityQuery,
    StructMotifQuery,
    StructureMotifResidue,
    Facet,
    Range as RCSBRange, # Renamed to avoid conflict with any built-in Range
    Attr
)
from rcsbsearchapi import rcsb_attributes as attrs
import json
import sys # Added import

# Attempt to import pypdb for the new tool
try:
    from pypdb import get_info as get_pdb_info_pypdb # Renamed to avoid conflict
    PYPDB_AVAILABLE = True
except ImportError:
    PYPDB_AVAILABLE = False
    print("Warning: pypdb not found. Protein information retrieval tool (get_protein_details_by_id_pypdb) will not be available.")
    # Define a dummy get_pdb_info_pypdb if pypdb is not available
    def get_pdb_info_pypdb(pdb_id: str, url_root: str = 'https://data.rcsb.org/rest/v1/co'): return None

# Create an MCP server
rcsb_host = os.getenv("RCSB_HOST", "0.0.0.0")
rcsb_port = int(os.getenv("RCSB_PORT", "8052"))
# Removed dependencies, host, and port from constructor
mcp = FastMCP("RCSB PDB MCP Server", description="Provides tools to query the RCSB PDB using rcsbsearchapi and pypdb for detailed protein information.")

# --- Text and Attribute Query Tools ---

@mcp.tool()
def text_search_pdb(query_string: str, return_type: str = "entry", results_verbosity: str = "compact", max_results: int = 10) -> str:
    """Perform a text search in RCSB PDB.
    Args:
        query_string: The text to search for (e.g., "hemoglobin").
        return_type: The type of identifiers to return (e.g., "entry", "polymer_entity", "assembly"). Default is "entry".
        results_verbosity: "compact" (IDs only), "minimal" (IDs and scores), "verbose" (all metadata). Default "compact".
        max_results: Maximum number of results to return. Default is 10.
    Returns:
        A JSON string with a list of entry IDs or dictionaries with metadata based on verbosity.
    
    Example Response Schema:
    for "compact" results_verbosity:
    {'data': ['2PGH',
        '3PEL',
        '3GOU',
        '6IHX',
        '1G08',
        '1G09',
        '1G0A',
        '2QSP',
        '5C6E',
        '3CIU']
    }

    for "minimal" results_verbosity:
    {'data': [{'identifier': '2PGH', 'score': 1.0},
            {'identifier': '3PEL', 'score': 0.9995626684955836},
            {'identifier': '3GOU', 'score': 0.9989625551654878},
            {'identifier': '6IHX', 'score': 0.9984230028463676},
            {'identifier': '1G08', 'score': 0.9978712997511598},
            {'identifier': '1G09', 'score': 0.9978712997511598},
            {'identifier': '1G0A', 'score': 0.9978712997511598},
            {'identifier': '2QSP', 'score': 0.9978712997511598},
            {'identifier': '5C6E', 'score': 0.9978712997511598},
            {'identifier': '3CIU', 'score': 0.9973749501914616}
    ]}

    for "verbose" results_verbosity:
    {'data': [{'identifier': '2PGH',
    'score': 1.0,
    'services': [{'service_type': 'full_text',
        'nodes': [{'node_id': 0,
        'original_score': 10.369708,
        'norm_score': 1.0}]}]},
    {'identifier': '3PEL',
    'score': 0.9995626684955836,
    'services': [{'service_type': 'full_text',
        'nodes': [{'node_id': 0,
        'original_score': 10.365173,
        'norm_score': 0.9995626684955836}]}]},
    {'identifier': '3GOU',
    'score': 0.9989625551654878,
    'services': [{'service_type': 'full_text',
        'nodes': [{'node_id': 0,
        'original_score': 10.35895,
        'norm_score': 0.9989625551654878}]}]},
    {'identifier': '6IHX',
    'score': 0.9984230028463676,
    'services': [{'service_type': 'full_text',
        'nodes': [{'node_id': 0,
        'original_score': 10.353355,
        'norm_score': 0.9984230028463676}]}]},
    {'identifier': '1G08',
    'score': 0.9978712997511598,
    'services': [{'service_type': 'full_text',
        'nodes': [{'node_id': 0,
        'original_score': 10.347634,
        'norm_score': 0.9978712997511598}]}]},
    {'identifier': '1G09',
    'score': 0.9978712997511598,
    'services': [{'service_type': 'full_text',
        'nodes': [{'node_id': 0,
        'original_score': 10.347634,
        'norm_score': 0.9978712997511598}]}]},
    {'identifier': '1G0A',
    'score': 0.9978712997511598,
    'services': [{'service_type': 'full_text',
        'nodes': [{'node_id': 0,
        'original_score': 10.347634,
        'norm_score': 0.9978712997511598}]}]},
    {'identifier': '2QSP',
    'score': 0.9978712997511598,
    'services': [{'service_type': 'full_text',
        'nodes': [{'node_id': 0,
        'original_score': 10.347634,
        'norm_score': 0.9978712997511598}]}]},
    {'identifier': '5C6E',
    'score': 0.9978712997511598,
    'services': [{'service_type': 'full_text',
        'nodes': [{'node_id': 0,
        'original_score': 10.347634,
        'norm_score': 0.9978712997511598}]}]},
    {'identifier': '3CIU',
    'score': 0.9973749501914616,
    'services': [{'service_type': 'full_text',
        'nodes': [{'node_id': 0,
        'original_score': 10.342487,
        'norm_score': 0.9973749501914616}]}]}]}
    """
    try:
        query = TextQuery(query_string)
        results = list(query(return_type=return_type, results_verbosity=results_verbosity))
        if len(results) > max_results:
            print(f"Query returned {len(results)} results. Truncating to max_results: {max_results}")
            results = results[:max_results]
        return json.dumps({"data": results})
    except Exception as e:
        return json.dumps({"error": f"Failed to perform text search for '{query_string}'.", "details": str(e)})

@mcp.tool()
def attribute_search_pdb(attribute_path: str, operator: str, value: list, return_type: str = "entry", results_verbosity: str = "compact", max_results: int = 10) -> str:
    """Perform an attribute search in RCSB PDB.
    Args:
        attribute_path: The path of the attribute to query (e.g., "exptl.method", "rcsb_entry_info.resolution_combined").
        operator: The comparison operator (e.g., "exact_match", "greater", "less_or_equal", "in").
        value: The value to compare against. For "in" operator, this should be a list.
        return_type: The type of identifiers to return. Default is "entry".
        results_verbosity: "compact", "minimal", or "verbose". Default "compact".
        max_results: Maximum number of results to return. Default is 10.
    Returns:
        A JSON string with a list of entry IDs or dictionaries with metadata.
    """
    try:
        query = AttributeQuery(attribute_path, operator, value)
        results = list(query(return_type=return_type, results_verbosity=results_verbosity))
        if len(results) > max_results:
            print(f"Query returned {len(results)} results. Truncating to max_results: {max_results}")
            results = results[:max_results]
        return json.dumps({"data": results})
    except Exception as e:
        return json.dumps({"error": f"Failed attribute search for {attribute_path} {operator} {value}.", "details": str(e)})

@mcp.tool()
def combined_text_and_attribute_search(text_query_string: str, attribute_filters: list[dict], logical_operator: str = "and", return_type: str = "entry", results_verbosity: str = "compact", max_results: int = 10) -> str:
    """Combines a text query with multiple attribute queries.
    Args:
        text_query_string: The main text query.
        attribute_filters: A list of attribute filter dictionaries. Each dict should have "attribute_path", "operator", and "value".
                           Example: [{"attribute_path": "rcsb_struct_symmetry.symbol", "operator": "exact_match", "value": "C2"}]
        logical_operator: How to combine the text query and attribute filters ("and" or "or"). Default "and".
        return_type: The type of identifiers to return. Default is "entry".
        results_verbosity: "compact", "minimal", or "verbose". Default "compact".
        max_results: Maximum number of results to return. Default is 10.
    Returns:
        A JSON string with a list of entry IDs or dictionaries with metadata.
    """
    try:
        main_query = TextQuery(text_query_string)
        
        attr_queries = []
        for filt in attribute_filters:
            attr_queries.append(AttributeQuery(filt["attribute_path"], filt["operator"], filt["value"]))

        if not attr_queries:
            final_query = main_query
        else:
            combined_attrs_query = attr_queries[0]
            if len(attr_queries) > 1:
                for i in range(1, len(attr_queries)):
                    combined_attrs_query = combined_attrs_query & attr_queries[i] # Default to AND for multiple attributes

            if logical_operator.lower() == "and":
                final_query = main_query & combined_attrs_query
            elif logical_operator.lower() == "or":
                final_query = main_query | combined_attrs_query
            else:
                return json.dumps({"error": "Invalid logical_operator. Must be 'and' or 'or'.", "details": f"Received: {logical_operator}"})
                
        results = list(final_query(return_type=return_type, results_verbosity=results_verbosity))
        if len(results) > max_results:
            print(f"Query returned {len(results)} results. Truncating to max_results: {max_results}")
            results = results[:max_results]
        return json.dumps({"data": results})
    except ValueError as ve:
        return json.dumps({"error": str(ve), "details": f"Logical operator was: {logical_operator}"})
    except Exception as e:
        return json.dumps({"error": "Failed combined text and attribute search.", "details": str(e)})

# --- Sequence Query Tools ---

@mcp.tool()
def sequence_identity_search(sequence: str, identity_cutoff: float = 0.9, e_value_cutoff: float = 1.0, sequence_type: str = "protein", return_type: str = "polymer_entity", max_results: int = 10) -> str:
    """Find PDB entities with sequence similarity.
    Args:
        sequence: The protein, DNA, or RNA sequence string.
        identity_cutoff: Minimum sequence identity (0.0 to 1.0). Default 0.9 (90%).
        e_value_cutoff: Maximum E-value for the match. Default 1.0.
        sequence_type: Type of sequence ("protein", "dna", "rna"). Default "protein".
        return_type: The type of identifiers to return. Default "polymer_entity".
        max_results: Maximum number of results to return. Default is 10.
    Returns:
        A JSON string with a list of polymer entity IDs.
    """
    try:
        query = SequenceQuery(sequence, sequence_type=sequence_type)
        results = list(query(return_type=return_type))
        if len(results) > max_results:
            print(f"Query returned {len(results)} results. Truncating to max_results: {max_results}")
            results = results[:max_results]
        return json.dumps({"data": results})
    except Exception as e:
        return json.dumps({"error": "Failed sequence identity search.", "details": str(e)})

@mcp.tool()
def sequence_motif_search(motif_pattern: str, pattern_type: str = "prosite", sequence_type: str = "protein", return_type: str = "polymer_entity", max_results: int = 10) -> str:
    """Search for sequences containing a specific motif.
    Args:
        motif_pattern: The motif pattern (e.g., "C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H." for PROSITE).
        pattern_type: Type of pattern ("prosite", "regex", "simple"). Default "prosite".
        sequence_type: Type of sequence ("protein", "dna", "rna"). Default "protein".
        return_type: The type of identifiers to return. Default "polymer_entity".
        max_results: Maximum number of results to return. Default is 10.
    Returns:
        A JSON string with a list of polymer entity IDs.
    """
    try:
        query = SeqMotifQuery(motif_pattern, pattern_type=pattern_type, sequence_type=sequence_type)
        results = list(query(return_type=return_type))
        if len(results) > max_results:
            print(f"Query returned {len(results)} results. Truncating to max_results: {max_results}")
            results = results[:max_results]
        return json.dumps({"data": results})
    except Exception as e:
        return json.dumps({"error": "Failed sequence motif search.", "details": str(e)})

# --- Structure Similarity Tools ---

@mcp.tool()
def structure_similarity_by_entry_id(entry_id: str, assembly_id: str = "1", operator: str = "strict_shape_match", target_search_space: str = "assembly", return_type: str = "assembly", max_results: int = 10) -> str:
    """Find structures similar to a given PDB entry.
    Args:
        entry_id: The PDB ID of the query structure (e.g., "4HHB").
        assembly_id: The assembly ID of the query structure. Default "1".
        operator: Similarity operator ("strict_shape_match" or "relaxed_shape_match"). Default "strict_shape_match".
        target_search_space: What to compare against ("assembly" or "polymer_entity_instance"). Default "assembly".
        return_type: The type of identifiers to return. Default "assembly".
        max_results: Maximum number of results to return. Default is 10.
    Returns:
        A JSON string with a list of assembly or polymer_entity_instance IDs.
    """
    try:
        query = StructSimilarityQuery(
            structure_search_type="entry_id",
            entry_id=entry_id,
            structure_input_type="assembly_id", # or chain_id if applicable
            assembly_id=assembly_id,
            operator=operator,
            target_search_space=target_search_space
        )
        results = list(query(return_type=return_type))
        if len(results) > max_results:
            print(f"Query returned {len(results)} results. Truncating to max_results: {max_results}")
            results = results[:max_results]
        return json.dumps({"data": results})
    except Exception as e:
        return json.dumps({"error": f"Failed structure similarity search for entry ID {entry_id}.", "details": str(e)})

@mcp.tool()
def structure_similarity_by_file_url(file_url: str, file_format: str, operator: str = "strict_shape_match", target_search_space: str = "assembly", return_type: str = "assembly", max_results: int = 10) -> str:
    """Find structures similar to a structure provided via a URL.
    Args:
        file_url: URL to the structure file (e.g., "https://files.rcsb.org/view/4HHB.cif").
        file_format: Format of the file ("cif", "bcif", "pdb", "cif.gz", "pdb.gz").
        operator: Similarity operator. Default "strict_shape_match".
        target_search_space: What to compare against. Default "assembly".
        return_type: The type of identifiers to return. Default "assembly".
        max_results: Maximum number of results to return. Default is 10.
    Returns:
        A JSON string with a list of assembly or polymer_entity_instance IDs.
    """
    try:
        query = StructSimilarityQuery(
            structure_search_type="file_url",
            file_url=file_url,
            file_format=file_format,
            operator=operator,
            target_search_space=target_search_space
        )
        results = list(query(return_type=return_type))
        if len(results) > max_results:
            print(f"Query returned {len(results)} results. Truncating to max_results: {max_results}")
            results = results[:max_results]
        return json.dumps({"data": results})
    except Exception as e:
        return json.dumps({"error": f"Failed structure similarity search for file URL {file_url}.", "details": str(e)})

# --- Structure Motif Tools ---
# Note: For file_upload, direct file path access from a server like this is complex and usually
# involves receiving the file content (e.g., as base64 string or multipart form data)
# and saving it temporarily. The rcsbsearchapi might expect a local file path.
# For simplicity, we'll skip direct file_upload from local machine path for now,
# as it's not straightforward in a typical MCP server context without more infrastructure.

# @mcp.tool()
# def structure_motif_search_by_entry_id(entry_id: str, residues: list[dict], backbone_distance_tolerance: int = 1, side_chain_distance_tolerance: int = 1, angle_tolerance: int = 1, rmsd_cutoff: float = 2.0, return_type: str = "polymer_entity", max_results: int = 10) -> str:
#     """Search for 3D structural motifs using a PDB entry as the reference for the motif.
#     Args:
#         entry_id: The PDB ID defining the motif (e.g., "2MNR").
#         residues: A list of residue definitions. Each dict must have "chain_id", "struct_oper_id" (operator, usually "1"),
#                   and "label_seq_id" (residue number). Optional "exchanges" (list of AA codes).
#                   Example: [{"chain_id": "A", "struct_oper_id": "1", "label_seq_id": 192, "exchanges": ["LYS", "HIS"]}]
#         backbone_distance_tolerance: Allowed backbone distance tolerance in Angstrom (0-3). Default 1.
#         side_chain_distance_tolerance: Allowed side-chain distance tolerance in Angstrom (0-3). Default 1.
#         angle_tolerance: Allowed angle tolerance in multiples of 20 degrees (0-3). Default 1.
#         rmsd_cutoff: Threshold above which hits will be filtered by RMSD (>=0). Default 2.0.
#         return_type: The type of identifiers to return. Default "polymer_entity".
#         max_results: Maximum number of results to return. Default is 10.
#     Returns:
#         A JSON string with a list of polymer entity IDs.
#     """
#     try:
#         motif_residues = [
#             StructureMotifResidue(
#                 chain_id=r.get("chain_id"),
#                 struct_oper_id=r.get("struct_oper_id", "1"), # Default to "1" if not provided
#                 label_seq_id=r.get("label_seq_id"),
#                 exchanges=r.get("exchanges")
#             ) for r in residues
#         ]
#         query = StructMotifQuery(
#             entry_id=entry_id,
#             residue_ids=motif_residues,
#             backbone_distance_tolerance=backbone_distance_tolerance,
#             side_chain_distance_tolerance=side_chain_distance_tolerance,
#             angle_tolerance=angle_tolerance,
#             rmsd_cutoff=rmsd_cutoff
#         )
#         results = list(query(return_type=return_type))
#         if len(results) > max_results:
#             print(f"Query returned {len(results)} results. Truncating to max_results: {max_results}")
#             results = results[:max_results]
#         return json.dumps({"data": results})
#     except Exception as e:
#         return json.dumps({"error": f"Failed structure motif search for entry ID {entry_id}.", "details": str(e)})

# --- PDB Protein Information Tool (using pypdb) ---

class PDBProteinInfoProvider:
    def __init__(self, pdb_id: str):
        self.pdb_id = pdb_id

    def _extract_simplified_pdb_data(self, pdb_data: dict) -> dict:
        """
        Extracts significant scientific properties from a complex PDB data dictionary.
        """
        if not pdb_data: # Handle case where pdb_data might be None (e.g. pypdb not found or invalid ID)
            return {"error": f"Could not retrieve data for PDB ID: '{self.pdb_id}'"}

        simplified_data = {}
        simplified_data['Pdb_Id'] = pdb_data.get('rcsb_id') or pdb_data.get('entry', {}).get('id')
        simplified_data['Title'] = pdb_data.get('struct', {}).get('title')
        authors_list = [author['name'] for author in pdb_data.get('audit_author', [])]
        simplified_data['Authors'] = ', '.join(authors_list) if authors_list else None
        citation = pdb_data.get('citation', [])
        if citation:
            primary_citation = next((c for c in citation if c.get('id') == 'primary'), citation[0])
            simplified_data['Journal'] = primary_citation.get('rcsb_journal_abbrev') or primary_citation.get('journal_abbrev')
            simplified_data['Year'] = primary_citation.get('year')
            simplified_data['Volume'] = primary_citation.get('journal_volume')
            page_first = primary_citation.get('page_first')
            page_last = primary_citation.get('page_last')
            if page_first and page_last:
                pages = f"{page_first}-{page_last}"
            elif page_first:
                pages = page_first
            else:
                pages = None
            simplified_data['Pages'] = pages
            simplified_data['Doi'] = primary_citation.get('pdbx_database_id_doi')
            simplified_data['Pubmed_Id'] = primary_citation.get('pdbx_database_id_pub_med')
        else:
            simplified_data['Journal'] = None
            simplified_data['Year'] = None
            simplified_data['Volume'] = None
            simplified_data['Pages'] = None
            simplified_data['Doi'] = None
            simplified_data['Pubmed_Id'] = None
        exptl = pdb_data.get('exptl', [])
        if exptl:
            simplified_data['Experiment_Method'] = exptl[0].get('method')
        else:
            simplified_data['Experiment_Method'] = None
        simplified_data['Molecular_Weight_(kDa)'] = pdb_data.get('rcsb_entry_info', {}).get('molecular_weight')
        simplified_data['Deposited_Model_Count'] = pdb_data.get('rcsb_entry_info', {}).get('deposited_model_count')
        simplified_data['Polymer_entity_count'] = pdb_data.get('rcsb_entry_info', {}).get('polymer_entity_count')
        simplified_data['Polymer_monomer_count'] = pdb_data.get('rcsb_entry_info', {}).get('deposited_polymer_monomer_count')
        simplified_data['Structural_Features'] = pdb_data.get('struct_keywords', {}).get('text')
        release_date = pdb_data.get('rcsb_accession_info', {}).get('initial_release_date', '')
        if 'T' in release_date:
            simplified_data['Release_Date'] = release_date.split('T')[0]
        else:
            simplified_data['Release_Date'] = release_date
        resolution = pdb_data.get('rcsb_entry_info', {}).get('resolution_combined', [None])
        if resolution and resolution[0] is not None:
            simplified_data['Resolution'] = resolution[0]
        else:
            simplified_data['Resolution'] = None
        return simplified_data

    def get_info(self) -> dict:
        if not PYPDB_AVAILABLE:
            return {"error": "pypdb is not installed. Cannot provide protein information."}
        try:
            raw_pdb_data = get_pdb_info_pypdb(self.pdb_id) # Use the renamed import
            if not raw_pdb_data: # Check if pypdb returned None (e.g. invalid PDB ID)
                 return {"error": f"No data found for PDB ID: '{self.pdb_id}'. It might be an invalid ID."}
            return self._extract_simplified_pdb_data(raw_pdb_data)
        except Exception as e:
            return {"error": f"Failed to retrieve or process PDB data for '{self.pdb_id}': {str(e)}"}

@mcp.tool()
def get_protein_details_by_id_pypdb(pdb_id_string: str) -> str:
    """
    Retrieves detailed information about a specific protein from its PDB ID using pypdb.
    This tool complements the rcsbsearchapi tools by providing a different set of details via pypdb.
    Args:
        pdb_id_string: The PDB ID of the protein (e.g., '6M0J', '1TIM').
    Returns:
        A JSON string containing the protein's details, or an error message.
    """
    if not PYPDB_AVAILABLE:
        return json.dumps({"error": "pypdb backend not available on the server for this tool."})
    
    provider = PDBProteinInfoProvider(pdb_id_string)
    info = provider.get_info()
    return json.dumps(info)

# --- Main execution for direct run or mcp dev ---
if __name__ == "__main__":
    print(f"Starting RCSB PDB MCP Server...")
    print(f"Server Name: {mcp.name}")
    sys.stdout.flush()

    # Default to SSE transport with host and port passed directly
    print(f"Attempting to run RCSB PDB MCP Server with FastMCP SSE transport on host {rcsb_host}, port {rcsb_port}")
    sys.stdout.flush()
    mcp.run(transport="sse", host=rcsb_host, port=rcsb_port)