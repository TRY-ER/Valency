import os
import sys 
from env_loader import load_env_vars 

load_env_vars() 

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
    Range as RCSBRange, 
    Attr
)
from rcsbsearchapi import rcsb_attributes as attrs
import json

# Create an MCP server
rcsb_host = os.getenv("RCSB_HOST", "0.0.0.0")
rcsb_port = int(os.getenv("RCSB_PORT", "8052"))
mcp = FastMCP(
    "RCSB PDB MCP Server", 
    description="Provides tools to query the RCSB PDB using rcsbsearchapi."
    # Removed dependencies argument
    # host and port settings are not used here, passed to run() for SSE
)

# --- Text and Attribute Query Tools ---

@mcp.tool()
def text_search_pdb(query_string: str, return_type: str = "entry", results_verbosity: str = "compact") -> str:
    """Perform a text search in RCSB PDB.
    Args:
        query_string: The text to search for (e.g., "hemoglobin").
        return_type: The type of identifiers to return (e.g., "entry", "polymer_entity", "assembly"). Default is "entry".
        results_verbosity: "compact" (IDs only), "minimal" (IDs and scores), "verbose" (all metadata). Default "compact".
    Returns:
        A JSON string with a list of entry IDs or dictionaries with metadata based on verbosity.
    """
    try:
        query = TextQuery(query_string)
        results = list(query(return_type=return_type, results_verbosity=results_verbosity))
        return json.dumps({"data": results})
    except Exception as e:
        return json.dumps({"error": f"Failed to perform text search for '{query_string}'.", "details": str(e)})

@mcp.tool()
def attribute_search_pdb(attribute_path: str, operator: str, value: list, return_type: str = "entry", results_verbosity: str = "compact") -> str:
    """Perform an attribute search in RCSB PDB.
    Args:
        attribute_path: The path of the attribute to query (e.g., "exptl.method", "rcsb_entry_info.resolution_combined").
        operator: The comparison operator (e.g., "exact_match", "greater", "less_or_equal", "in").
        value: The value to compare against. For "in" operator, this should be a list.
        return_type: The type of identifiers to return. Default is "entry".
        results_verbosity: "compact", "minimal", or "verbose". Default "compact".
    Returns:
        A JSON string with a list of entry IDs or dictionaries with metadata.
    """
    try:
        query = AttributeQuery(attribute_path, operator, value)
        results = list(query(return_type=return_type, results_verbosity=results_verbosity))
        return json.dumps({"data": results})
    except Exception as e:
        return json.dumps({"error": f"Failed attribute search for {attribute_path} {operator} {value}.", "details": str(e)})

@mcp.tool()
def combined_text_and_attribute_search(text_query_string: str, attribute_filters: list[dict], logical_operator: str = "and", return_type: str = "entry", results_verbosity: str = "compact") -> str:
    """Combines a text query with multiple attribute queries.
    Args:
        text_query_string: The main text query.
        attribute_filters: A list of attribute filter dictionaries. Each dict should have "attribute_path", "operator", and "value".
                           Example: [{"attribute_path": "rcsb_struct_symmetry.symbol", "operator": "exact_match", "value": "C2"}]
        logical_operator: How to combine the text query and attribute filters ("and" or "or"). Default "and".
        return_type: The type of identifiers to return. Default is "entry".
        results_verbosity: "compact", "minimal", or "verbose". Default "compact".
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
        return json.dumps({"data": results})
    except ValueError as ve:
        return json.dumps({"error": str(ve), "details": f"Logical operator was: {logical_operator}"})
    except Exception as e:
        return json.dumps({"error": "Failed combined text and attribute search.", "details": str(e)})

# --- Sequence Query Tools ---

@mcp.tool()
def sequence_identity_search(sequence: str, identity_cutoff: float = 0.9, e_value_cutoff: float = 1.0, sequence_type: str = "protein", return_type: str = "polymer_entity") -> str:
    """Find PDB entities with sequence similarity.
    Args:
        sequence: The protein, DNA, or RNA sequence string.
        identity_cutoff: Minimum sequence identity (0.0 to 1.0). Default 0.9 (90%).
        e_value_cutoff: Maximum E-value for the match. Default 1.0.
        sequence_type: Type of sequence ("protein", "dna", "rna"). Default "protein".
        return_type: The type of identifiers to return. Default "polymer_entity".
    Returns:
        A JSON string with a list of polymer entity IDs.
    """
    try:
        query = SequenceQuery(sequence, e_value_cutoff=e_value_cutoff, identity_cutoff=identity_cutoff, sequence_type=sequence_type)
        results = list(query(return_type=return_type))
        return json.dumps({"data": results})
    except Exception as e:
        return json.dumps({"error": "Failed sequence identity search.", "details": str(e)})

@mcp.tool()
def sequence_motif_search(motif_pattern: str, pattern_type: str = "prosite", sequence_type: str = "protein", return_type: str = "polymer_entity") -> str:
    """Search for sequences containing a specific motif.
    Args:
        motif_pattern: The motif pattern (e.g., "C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H." for PROSITE).
        pattern_type: Type of pattern ("prosite", "regex", "simple"). Default "prosite".
        sequence_type: Type of sequence ("protein", "dna", "rna"). Default "protein".
        return_type: The type of identifiers to return. Default "polymer_entity".
    Returns:
        A JSON string with a list of polymer entity IDs.
    """
    try:
        query = SeqMotifQuery(motif_pattern, pattern_type=pattern_type, sequence_type=sequence_type)
        results = list(query(return_type=return_type))
        return json.dumps({"data": results})
    except Exception as e:
        return json.dumps({"error": "Failed sequence motif search.", "details": str(e)})

# --- Structure Similarity Tools ---

@mcp.tool()
def structure_similarity_by_entry_id(entry_id: str, assembly_id: str = "1", operator: str = "strict_shape_match", target_search_space: str = "assembly", return_type: str = "assembly") -> str:
    """Find structures similar to a given PDB entry.
    Args:
        entry_id: The PDB ID of the query structure (e.g., "4HHB").
        assembly_id: The assembly ID of the query structure. Default "1".
        operator: Similarity operator ("strict_shape_match" or "relaxed_shape_match"). Default "strict_shape_match".
        target_search_space: What to compare against ("assembly" or "polymer_entity_instance"). Default "assembly".
        return_type: The type of identifiers to return. Default "assembly".
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
        return json.dumps({"data": results})
    except Exception as e:
        return json.dumps({"error": f"Failed structure similarity search for entry ID {entry_id}.", "details": str(e)})

@mcp.tool()
def structure_similarity_by_file_url(file_url: str, file_format: str, operator: str = "strict_shape_match", target_search_space: str = "assembly", return_type: str = "assembly") -> str:
    """Find structures similar to a structure provided via a URL.
    Args:
        file_url: URL to the structure file (e.g., "https://files.rcsb.org/view/4HHB.cif").
        file_format: Format of the file ("cif", "bcif", "pdb", "cif.gz", "pdb.gz").
        operator: Similarity operator. Default "strict_shape_match".
        target_search_space: What to compare against. Default "assembly".
        return_type: The type of identifiers to return. Default "assembly".
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
        return json.dumps({"data": results})
    except Exception as e:
        return json.dumps({"error": f"Failed structure similarity search for file URL {file_url}.", "details": str(e)})

# --- Structure Motif Tools ---
# Note: For file_upload, direct file path access from a server like this is complex and usually
# involves receiving the file content (e.g., as base64 string or multipart form data)
# and saving it temporarily. The rcsbsearchapi might expect a local file path.
# For simplicity, we'll skip direct file_upload from local machine path for now,
# as it's not straightforward in a typical MCP server context without more infrastructure.

@mcp.tool()
def structure_motif_search_by_entry_id(entry_id: str, residues: list[dict], backbone_distance_tolerance: int = 1, side_chain_distance_tolerance: int = 1, angle_tolerance: int = 1, rmsd_cutoff: float = 2.0, return_type: str = "polymer_entity") -> str:
    """Search for 3D structural motifs using a PDB entry as the reference for the motif.
    Args:
        entry_id: The PDB ID defining the motif (e.g., "2MNR").
        residues: A list of residue definitions. Each dict must have "chain_id", "struct_oper_id" (operator, usually "1"),
                  and "label_seq_id" (residue number). Optional "exchanges" (list of AA codes).
                  Example: [{"chain_id": "A", "struct_oper_id": "1", "label_seq_id": 192, "exchanges": ["LYS", "HIS"]}]
        backbone_distance_tolerance: Allowed backbone distance tolerance in Angstrom (0-3). Default 1.
        side_chain_distance_tolerance: Allowed side-chain distance tolerance in Angstrom (0-3). Default 1.
        angle_tolerance: Allowed angle tolerance in multiples of 20 degrees (0-3). Default 1.
        rmsd_cutoff: Threshold above which hits will be filtered by RMSD (>=0). Default 2.0.
        return_type: The type of identifiers to return. Default "polymer_entity".
    Returns:
        A JSON string with a list of polymer entity IDs.
    """
    try:
        motif_residues = [
            StructureMotifResidue(
                chain_id=r.get("chain_id"),
                struct_oper_id=r.get("struct_oper_id", "1"), # Default to "1" if not provided
                label_seq_id=r.get("label_seq_id"),
                exchanges=r.get("exchanges")
            ) for r in residues
        ]
        query = StructMotifQuery(
            entry_id=entry_id,
            residue_ids=motif_residues,
            backbone_distance_tolerance=backbone_distance_tolerance,
            side_chain_distance_tolerance=side_chain_distance_tolerance,
            angle_tolerance=angle_tolerance,
            rmsd_cutoff=rmsd_cutoff
        )
        results = list(query(return_type=return_type))
        return json.dumps({"data": results})
    except Exception as e:
        return json.dumps({"error": f"Failed structure motif search for entry ID {entry_id}.", "details": str(e)})

# --- Facet/Aggregation Tools ---

@mcp.tool()
def get_term_facets(attribute_query_dict: dict, facet_name: str, facet_attribute: str, min_interval_population: int = 1, max_num_intervals: int = 10) -> str:
    """Perform a query and get term-based aggregations (facets).
    Args:
        attribute_query_dict: Dictionary defining the base attribute query (e.g., {"attribute_path": "rcsb_entry_info.structure_determination_methodology", "operator": "exact_match", "value": "experimental"}).
        facet_name: A descriptive name for this facet (e.g., "Experimental Methods").
        facet_attribute: The attribute to aggregate on (e.g., "exptl.method").
        min_interval_population: Minimum count for a bucket to be returned. Default 1.
        max_num_intervals: Maximum number of buckets to return. Default 10.
    Returns:
        A JSON string with a list of dictionaries, each representing a facet bucket.
    """
    try:
        base_query = AttributeQuery(attribute_query_dict["attribute_path"], attribute_query_dict["operator"], attribute_query_dict["value"])
        facet_obj = Facet(name=facet_name, aggregation_type="terms", attribute=facet_attribute, min_interval_population=min_interval_population, max_num_intervals=max_num_intervals)
        
        # The .facets attribute is accessed after calling the query object with the facet
        results = base_query(facets=facet_obj).facets
        return json.dumps({"data": results})
    except Exception as e:
        return json.dumps({"error": "Failed to get term facets.", "details": str(e)})


@mcp.tool()
def get_histogram_facets(attribute_query_dict: dict, facet_name: str, facet_attribute: str, interval: float, min_interval_population: int = 1, return_type: str = "entry") -> str:
    """Perform a query and get histogram-based aggregations for numeric attributes.
    Args:
        attribute_query_dict: Dictionary for the base attribute query.
        facet_name: Descriptive name for the facet (e.g., "Formula Weight Distribution").
        facet_attribute: The numeric attribute for histogram (e.g., "rcsb_polymer_entity.formula_weight").
        interval: The size of the histogram intervals.
        min_interval_population: Minimum count for a bucket. Default 1.
        return_type: Return type for the base query. Default "entry".
    Returns:
        A JSON string with a list of dictionaries, each representing a histogram bucket.
    """
    try:
        base_query = AttributeQuery(attribute_query_dict["attribute_path"], attribute_query_dict["operator"], attribute_query_dict["value"])
        facet_obj = Facet(name=facet_name, aggregation_type="histogram", attribute=facet_attribute, interval=interval, min_interval_population=min_interval_population)
        results = base_query(return_type=return_type, facets=facet_obj).facets
        return json.dumps({"data": results})
    except Exception as e:
        return json.dumps({"error": "Failed to get histogram facets.", "details": str(e)})

@mcp.tool()
def get_date_histogram_facets(attribute_query_dict: dict, facet_name: str, facet_attribute: str, interval: str = "year", min_interval_population: int = 1) -> str:
    """Perform a query and get date histogram aggregations.
    Args:
        attribute_query_dict: Dictionary for the base attribute query.
        facet_name: Descriptive name for the facet (e.g., "Release Dates by Year").
        facet_attribute: The date attribute for histogram (e.g., "rcsb_accession_info.initial_release_date").
        interval: Interval for date histogram (must be "year" as per notebook, though API might support others). Default "year".
        min_interval_population: Minimum count for a bucket. Default 1.
    Returns:
        A JSON string with a list of dictionaries, each representing a date histogram bucket.
    """
    try:
        base_query = AttributeQuery(attribute_query_dict["attribute_path"], attribute_query_dict["operator"], attribute_query_dict["value"])
        facet_obj = Facet(name=facet_name, aggregation_type="date_histogram", attribute=facet_attribute, interval=interval, min_interval_population=min_interval_population)
        results = base_query(facets=facet_obj).facets
        return json.dumps({"data": results})
    except Exception as e:
        return json.dumps({"error": "Failed to get date histogram facets.", "details": str(e)})

@mcp.tool()
def get_cardinality_facet(attribute_query_dict: dict, facet_name: str, facet_attribute: str, precision_threshold: int = 40000) -> str:
    """Get the count of distinct values for a field.
    Args:
        attribute_query_dict: Dictionary for the base attribute query.
        facet_name: Descriptive name for the facet (e.g., "Unique Organism Count").
        facet_attribute: The attribute for which to count distinct values (e.g., "rcsb_entity_source_organism.ncbi_scientific_name").
        precision_threshold: Precision threshold for cardinality calculation. Default 40000.
    Returns:
        A JSON string containing a single dictionary with the cardinality result.
    """
    try:
        base_query = AttributeQuery(attribute_query_dict["attribute_path"], attribute_query_dict["operator"], attribute_query_dict["value"])
        facet_obj = Facet(name=facet_name, aggregation_type="cardinality", attribute=facet_attribute, precision_threshold=precision_threshold)
        results = base_query(facets=facet_obj).facets
        return json.dumps({"data": results})
    except Exception as e:
        return json.dumps({"error": "Failed to get cardinality facet.", "details": str(e)})
    
# --- Utility Tools ---

@mcp.tool()
def count_query_results(query_type: str, query_params: dict, return_content_type: list[str] | None = None) -> str:
    """Counts the number of results for a given query.
    Args:
        query_type: Type of query ("text", "attribute", "sequence", "seq_motif", "struct_similarity", "struct_motif").
        query_params: Dictionary of parameters for the specified query type.
                       For "text": {"query_string": "value"}
                       For "attribute": {"attribute_path": "path", "operator": "op", "value": "val"}
                       ... and so on, matching the primary args of other tools.
        return_content_type: Optional list to specify content types like ["computational", "experimental"].
    Returns:
        A JSON string with the number of results.
    """
    try:
        query = None
        if query_type == "text":
            query = TextQuery(query_params["query_string"])
        elif query_type == "attribute":
            query = AttributeQuery(query_params["attribute_path"], query_params["operator"], query_params["value"])
        # Add other query types if needed, ensuring params match their constructors
        else:
            return json.dumps({"error": f"Unsupported query_type for counting: {query_type}"})

        if return_content_type:
            count = query(return_counts=True, return_content_type=return_content_type)
        else:
            count = query(return_counts=True)
        return json.dumps({"count": count})
    except KeyError as ke:
        return json.dumps({"error": f"Missing parameter for query_type '{query_type}'.", "details": str(ke)})
    except ValueError as ve:
        return json.dumps({"error": str(ve), "details": f"Query type was: {query_type}"})
    except Exception as e:
        return json.dumps({"error": "Failed to count query results.", "details": str(e)})


# --- Main execution for direct run or mcp dev ---
if __name__ == "__main__":
    # transport = os.getenv("MCP_TRANSPORT", "sse") # No longer needed for direct SSE run
    print(f"Starting RCSB PDB MCP Server...")
    print(f"Server Name: {mcp.name}")
    sys.stdout.flush() 

    # Default to SSE transport with host and port passed directly
    print(f"Attempting to run RCSB PDB MCP Server with FastMCP SSE transport on host {rcsb_host}, port {rcsb_port}")
    sys.stdout.flush() 
    mcp.run(transport="sse", host=rcsb_host, port=rcsb_port)
    # else:
    #     raise ValueError(f"Unknown transport: {transport}. Supported: 'stdio', 'sse'.")

