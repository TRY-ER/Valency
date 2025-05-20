import os
from env_loader import load_env_vars # Changed from dotenv to env_loader

load_env_vars() # Load variables from .env file using custom loader

from mcp.server.fastmcp import FastMCP
from chembl_webresource_client.new_client import new_client
from chembl_webresource_client.utils import utils as chembl_utils
import json

# Create an MCP server
chembl_host = os.getenv("CHEMBL_HOST", "0.0.0.0")
chembl_port = int(os.getenv("CHEMBL_PORT", "8051"))
mcp = FastMCP("ChEMBL MCP Server", dependencies=["chembl_webresource_client"], host=chembl_host, port=chembl_port)

# --- Molecule Tools ---

@mcp.tool()
def get_molecule_by_chembl_id(chembl_id: str, only_fields: list[str] | None = None) -> str:
    """Retrieve a specific molecule by its unique ChEMBL ID (e.g., 'CHEMBL192').
    This is the most direct way to get a molecule if its ChEMBL ID is known.
    Args:
        chembl_id: The ChEMBL ID of the molecule to retrieve.
        only_fields: Optional list of specific fields to return (e.g., ['molecule_chembl_id', 'pref_name', 'molecule_structures']).
                     If None, all available fields are returned.
    Returns:
        A JSON string containing the molecule data, or None if not found.
    """
    try:
        molecule_client = new_client.molecule
        if only_fields:
            results = molecule_client.filter(chembl_id=chembl_id).only(only_fields)
        else:
            results = molecule_client.filter(chembl_id=chembl_id)
        return json.dumps({"data": results[0] if results else None})
    except Exception as e:
        return json.dumps({"error": f"Failed to get molecule by ChEMBL ID '{chembl_id}'.", "details": str(e)})

@mcp.tool()
def find_molecule_by_pref_name(pref_name: str, exact_match: bool = True, only_fields: list[str] | None = None) -> str:
    """Search for molecules by their preferred name (e.g., 'Aspirin').
    Useful when the exact ChEMBL ID is not known but the common name is.
    Args:
        pref_name: The preferred name to search for.
        exact_match: If True (default), performs an exact case-insensitive match. 
                     If False, performs a case-insensitive partial match (contains).
        only_fields: Optional list of specific fields to return.
    Returns:
        A JSON string containing a list of dictionaries, each representing a matching molecule.
    """
    try:
        molecule_client = new_client.molecule
        if exact_match:
            query = molecule_client.filter(pref_name__iexact=pref_name)
        else:
            query = molecule_client.filter(pref_name__icontains=pref_name)
        if only_fields:
            query = query.only(only_fields)
        return json.dumps({"data": list(query)})
    except Exception as e:
        return json.dumps({"error": f"Failed to find molecule by preferred name '{pref_name}'.", "details": str(e)})

@mcp.tool()
def find_molecule_by_synonym(synonym: str, only_fields: list[str] | None = None) -> str:
    """Find molecules by their synonyms (e.g., 'Viagra' for Sildenafil).
    Helpful if a molecule is known by an alternative name or brand name.
    Args:
        synonym: The synonym to search for (case-insensitive exact match).
        only_fields: Optional list of specific fields to return.
    Returns:
        A JSON string containing a list of dictionaries, each representing a matching molecule.
    """
    try:
        molecule_client = new_client.molecule
        query = molecule_client.filter(molecule_synonyms__molecule_synonym__iexact=synonym)
        if only_fields:
            query = query.only(only_fields)
        return json.dumps({"data": list(query)})
    except Exception as e:
        return json.dumps({"error": f"Failed to find molecule by synonym '{synonym}'.", "details": str(e)})

@mcp.tool()
def get_molecules_by_chembl_ids(chembl_ids: list[str], only_fields: list[str] | None = None) -> str:
    """Retrieve multiple molecules by providing a list of their ChEMBL IDs.
    Efficient for batch lookups.
    Args:
        chembl_ids: A list of ChEMBL IDs (e.g., ['CHEMBL25', 'CHEMBL192']).
        only_fields: Optional list of specific fields to return for each molecule.
    Returns:
        A JSON string containing a list of dictionaries, each representing a molecule.
    """
    try:
        molecule_client = new_client.molecule
        query = molecule_client.filter(molecule_chembl_id__in=chembl_ids)
        if only_fields:
            query = query.only(only_fields)
        return json.dumps({"data": list(query)})
    except Exception as e:
        return json.dumps({"error": f"Failed to get molecules by ChEMBL IDs: {chembl_ids}.", "details": str(e)})

@mcp.tool()
def get_molecule_image_svg(chembl_id: str) -> str:
    """Get the 2D chemical structure image of a molecule as an SVG string.
    Useful for visualizing a molecule when its ChEMBL ID is known.
    Args:
        chembl_id: The ChEMBL ID of the molecule.
    Returns:
        A JSON string containing an SVG string representing the molecule's image, or None if not found/available.
    """
    try:
        image_client = new_client.image
        image_client.set_format('svg')
        svg_data = image_client.get(chembl_id)
        return json.dumps({"svg_image": svg_data})
    except Exception as e:
        return json.dumps({"error": f"Failed to get molecule image for ChEMBL ID '{chembl_id}'.", "details": str(e)})

@mcp.tool()
def find_similar_molecules_by_smiles(smiles: str, similarity_threshold: int = 70, only_fields: list[str] | None = None) -> str:
    """Find molecules structurally similar to a given molecule represented by its SMILES string.
    Uses Tanimoto similarity. A threshold of 70 means 70% similar or more.
    Args:
        smiles: The SMILES string of the query molecule (e.g., 'CCO').
        similarity_threshold: The minimum similarity percentage (0-100, default 70).
        only_fields: Optional list of specific fields to return for similar molecules.
    Returns:
        A JSON string containing a list of dictionaries, each representing a similar molecule, including its similarity score.
    """
    try:
        similarity_client = new_client.similarity
        query = similarity_client.filter(smiles=smiles, similarity=similarity_threshold)
        if only_fields:
            query = query.only(only_fields)
        return json.dumps({"data": list(query)})
    except Exception as e:
        return json.dumps({"error": f"Failed to find similar molecules by SMILES '{smiles}'.", "details": str(e)})

@mcp.tool()
def find_similar_molecules_by_chembl_id(chembl_id: str, similarity_threshold: int = 70, only_fields: list[str] | None = None) -> str:
    """Find molecules structurally similar to a given molecule identified by its ChEMBL ID.
    Uses Tanimoto similarity. A threshold of 70 means 70% similar or more.
    Args:
        chembl_id: The ChEMBL ID of the query molecule.
        similarity_threshold: The minimum similarity percentage (0-100, default 70).
        only_fields: Optional list of specific fields to return for similar molecules.
    Returns:
        A JSON string containing a list of dictionaries, each representing a similar molecule, including its similarity score.
    """
    try:
        similarity_client = new_client.similarity
        query = similarity_client.filter(chembl_id=chembl_id, similarity=similarity_threshold)
        if only_fields:
            query = query.only(only_fields)
        return json.dumps({"data": list(query)})
    except Exception as e:
        return json.dumps({"error": f"Failed to find similar molecules by ChEMBL ID '{chembl_id}'.", "details": str(e)})

@mcp.tool()
def get_approved_drugs(order_by_mw: bool = False, only_fields: list[str] | None = None) -> str:
    """Retrieve all drugs that have reached the maximum clinical trial phase (approved drugs).
    Args:
        order_by_mw: If True, sorts the results by molecular weight (ascending).
        only_fields: Optional list of specific fields to return.
    Returns:
        A JSON string containing a list of dictionaries, each representing an approved drug.
    """
    try:
        molecule_client = new_client.molecule
        query = molecule_client.filter(max_phase=4)
        if order_by_mw:
            query = query.order_by('molecule_properties__mw_freebase')
        if only_fields:
            query = query.only(only_fields)
        return json.dumps({"data": list(query)})
    except Exception as e:
        return json.dumps({"error": "Failed to get approved drugs.", "details": str(e)})

# --- Activity Tools ---

@mcp.tool()
def get_activities_for_target(target_chembl_id: str, standard_type: str | None = "IC50", only_fields: list[str] | None = None) -> str:
    """Fetch bioactivity data (e.g., IC50, Ki) for a specific biological target.
    Useful for finding out which compounds are active against a particular target.
    Args:
        target_chembl_id: The ChEMBL ID of the target (e.g., 'CHEMBL240' for EGFR).
        standard_type: The type of bioactivity measurement (e.g., 'IC50', 'Ki', 'EC50'). Default is 'IC50'. 
                       Pass None to get all activity types.
        only_fields: Optional list of specific fields to return for each activity record.
    Returns:
        A JSON string containing a list of dictionaries, each representing an activity record.
    """
    try:
        activity_client = new_client.activity
        query = activity_client.filter(target_chembl_id=target_chembl_id)
        if standard_type:
            query = query.filter(standard_type__iexact=standard_type)
        if only_fields:
            query = query.only(only_fields)
        return json.dumps({"data": list(query)})
    except Exception as e:
        return json.dumps({"error": f"Failed to get activities for target ChEMBL ID '{target_chembl_id}'.", "details": str(e)})

@mcp.tool()
def get_activities_for_molecule(molecule_chembl_id: str, pchembl_value_exists: bool = True, only_fields: list[str] | None = None) -> str:
    """Retrieve all recorded bioactivities for a specific molecule.
    Useful for understanding the biological profile of a compound.
    Args:
        molecule_chembl_id: The ChEMBL ID of the molecule.
        pchembl_value_exists: If True (default), only returns activities that have a pChEMBL value (a standardized measure of potency).
        only_fields: Optional list of specific fields to return.
    Returns:
        A JSON string containing a list of dictionaries, each representing an activity record.
    """
    try:
        activity_client = new_client.activity
        query = activity_client.filter(molecule_chembl_id=molecule_chembl_id)
        if pchembl_value_exists:
            query = query.filter(pchembl_value__isnull=False)
        if only_fields:
            query = query.only(only_fields)
        return json.dumps({"data": list(query)})
    except Exception as e:
        return json.dumps({"error": f"Failed to get activities for molecule ChEMBL ID '{molecule_chembl_id}'.", "details": str(e)})

# --- Target Tools ---

@mcp.tool()
def find_target_by_gene_name(gene_name: str, only_fields: list[str] | None = None) -> str:
    """Search for biological targets (e.g., proteins, protein families) by a gene name or symbol.
    This tool searches within target synonyms, so it can find targets even if the gene name is not the preferred name.
    Args:
        gene_name: The gene name or symbol to search for (e.g., 'BRCA1', 'EGFR'). Case-insensitive contains match.
        only_fields: Optional list of specific fields to return (e.g., ['organism', 'pref_name', 'target_type', 'target_chembl_id']).
    Returns:
        A JSON string containing a list of dictionaries, each representing a matching target.
    """
    try:
        target_client = new_client.target
        query = target_client.filter(target_synonym__icontains=gene_name)
        if only_fields:
            query = query.only(only_fields)
        return json.dumps({"data": list(query)})
    except Exception as e:
        return json.dumps({"error": f"Failed to find target by gene name '{gene_name}'.", "details": str(e)})

# --- Generic Filter Tools ---

@mcp.tool()
def get_molecules_by_filter(filters: dict[str, any], only_fields: list[str] | None = None, order_by: list[str] | None = None) -> str:
    """A flexible tool to retrieve molecules based on a custom set of filter conditions.
    Use this for complex queries not covered by more specific tools.
    Args:
        filters: A dictionary where keys are ChEMBL molecule fields (e.g., 'molecule_properties__mw_freebase') 
                 and values are the conditions (e.g., {'molecule_properties__mw_freebase__gte': 200, 'max_phase': 4}).
                 Supports Django-style lookups (e.g., '__gte', '__icontains').
        only_fields: Optional list of specific fields to return.
        order_by: Optional list of fields to sort the results by (e.g., ['molecule_properties__mw_freebase']). 
                  Prefix with '-' for descending order (e.g., ['-num_ro5_violations']).
    Returns:
        A JSON string containing a list of dictionaries, each representing a molecule matching the criteria.
    Example Usage:
        To get molecules with molecular weight >= 200, ALOGP <= 5, and max_phase = 4:
        get_molecules_by_filter(filters={"molecule_properties__mw_freebase__gte": 200, 
                                       "molecule_properties__alogp__lte": 5, 
                                       "max_phase": 4})
    """
    try:
        molecule_client = new_client.molecule
        query = molecule_client.filter(**filters)
        if order_by:
            query = query.order_by(*order_by)
        if only_fields:
            query = query.only(only_fields)
        return json.dumps({"data": list(query)})
    except Exception as e:
        return json.dumps({"error": f"Failed to get molecules by filter: {filters}.", "details": str(e)})

@mcp.tool()
def get_activities_by_filter(filters: dict[str, any], only_fields: list[str] | None = None, order_by: list[str] | None = None) -> str:
    """A flexible tool to retrieve bioactivity data based on a custom set of filter conditions.
    Use this for complex queries on activity data.
    Args:
        filters: A dictionary of filter conditions for activity fields (e.g., {'pchembl_value__gte': 6.0, 'standard_type__iexact': 'IC50'}).
        only_fields: Optional list of specific fields to return.
        order_by: Optional list of fields to sort the results by (e.g., ['-pchembl_value']).
    Returns:
        A JSON string containing a list of dictionaries, each representing an activity record.
    Example Usage:
        To get IC50 activities with pChEMBL value >= 6.0 for target CHEMBL240:
        get_activities_by_filter(filters={"pchembl_value__gte": 6.0, 
                                          "standard_type__iexact": "IC50", 
                                          "target_chembl_id": "CHEMBL240"})
    """
    try:
        activity_client = new_client.activity
        query = activity_client.filter(**filters)
        if order_by:
            query = query.order_by(*order_by)
        if only_fields:
            query = query.only(only_fields)
        return json.dumps({"data": list(query)})
    except Exception as e:
        return json.dumps({"error": f"Failed to get activities by filter: {filters}.", "details": str(e)})

@mcp.tool()
def get_targets_by_filter(filters: dict[str, any], only_fields: list[str] | None = None, order_by: list[str] | None = None) -> str:
    """A flexible tool to retrieve biological targets based on a custom set of filter conditions.
    Use this for complex queries on target data.
    Args:
        filters: A dictionary of filter conditions for target fields (e.g., {'pref_name__startswith': 'Cytochrome P450', 'target_type': 'SINGLE PROTEIN'}).
        only_fields: Optional list of specific fields to return.
        order_by: Optional list of fields to sort the results by (e.g., ['pref_name']).
    Returns:
        A JSON string containing a list of dictionaries, each representing a target matching the criteria.
    Example Usage:
        To get human single protein targets whose preferred name starts with 'Cytochrome P450':
        get_targets_by_filter(filters={"pref_name__startswith": "Cytochrome P450", 
                                       "target_type": "SINGLE PROTEIN", 
                                       "organism__istartswith": "Homo sapiens"})
    """
    try:
        target_client = new_client.target
        query = target_client.filter(**filters)
        if order_by:
            query = query.order_by(*order_by)
        if only_fields:
            query = query.only(only_fields)
        return json.dumps({"data": list(query)})
    except Exception as e:
        return json.dumps({"error": f"Failed to get targets by filter: {filters}.", "details": str(e)})

# --- Utility Tools ---

@mcp.tool()
def smiles_to_ctab(smiles: str) -> str:
    """Convert a SMILES (Simplified Molecular Input Line Entry System) string to a CTAB (Chemical Table) block.
    CTAB is a text-based format for representing chemical structures.
    Args:
        smiles: The SMILES string to convert.
    Returns:
        A JSON string containing the molecule in CTAB format.
    """
    try:
        ctab_data = chembl_utils.smiles2ctab(smiles)
        return json.dumps({"ctab": ctab_data})
    except Exception as e:
        return json.dumps({"error": f"Failed to convert SMILES '{smiles}' to CTAB.", "details": str(e)})

@mcp.tool()
def compute_molecular_descriptors(smiles: str) -> str:
    """Calculate a set of physicochemical properties and descriptors for a molecule from its SMILES string.
    Examples of descriptors include Molecular Weight, AlogP, Number of Hydrogen Bond Donors/Acceptors, etc.
    Args:
        smiles: The SMILES string of the molecule.
    Returns:
        A JSON string containing various calculated molecular descriptors.
    """
    try:
        ctab = chembl_utils.smiles2ctab(smiles)
        descriptors_json_str = chembl_utils.chemblDescriptors(ctab)
        descriptors_list = json.loads(descriptors_json_str)
        return json.dumps({"data": descriptors_list[0] if descriptors_list else {}})
    except Exception as e:
        return json.dumps({"error": f"Failed to compute molecular descriptors for SMILES '{smiles}'.", "details": str(e)})

@mcp.tool()
def compute_structural_alerts(smiles: str) -> str:
    """Identify known structural alerts (e.g., toxicophores, groups associated with reactivity) within a molecule from its SMILES string.
    Useful for early-stage assessment of potential liabilities.
    Args:
        smiles: The SMILES string of the molecule.
    Returns:
        A JSON string containing a list of dictionaries, where each dictionary represents a structural alert found in the molecule.
    """
    try:
        ctab = chembl_utils.smiles2ctab(smiles)
        alerts_json_str = chembl_utils.structuralAlerts(ctab)
        alerts_list = json.loads(alerts_json_str)
        return json.dumps({"data": alerts_list[0] if alerts_list else []})
    except Exception as e:
        return json.dumps({"error": f"Failed to compute structural alerts for SMILES '{smiles}'.", "details": str(e)})

@mcp.tool()
def standardize_molecule_from_smiles(smiles: str) -> str:
    """Standardize a molecular structure provided as a SMILES string.
    This process typically involves neutralizing charges, removing salts, and applying a consistent representation.
    Args:
        smiles: The SMILES string of the molecule to standardize.
    Returns:
        A JSON string containing the standardized molecular information, usually including the standardized SMILES.
    """
    try:
        ctab = chembl_utils.smiles2ctab(smiles)
        standardized_json_str = chembl_utils.standardize(ctab)
        return json.dumps({"data": json.loads(standardized_json_str)})
    except Exception as e:
        return json.dumps({"error": f"Failed to standardize molecule from SMILES '{smiles}'.", "details": str(e)})

@mcp.tool()
def get_parent_molecule_from_smiles(smiles: str) -> str:
    """For a molecule that might be a salt or a mixture (represented by its SMILES string), this tool attempts to identify and return its parent structure.
    This often means removing counter-ions or selecting the largest covalent component.
    Args:
        smiles: The SMILES string of the molecule.
    Returns:
        A JSON string containing information about the parent molecule, usually including its SMILES.
    """
    try:
        ctab = chembl_utils.smiles2ctab(smiles)
        parent_json_str = chembl_utils.getParent(ctab)
        return json.dumps({"data": json.loads(parent_json_str)})
    except Exception as e:
        return json.dumps({"error": f"Failed to get parent molecule from SMILES '{smiles}'.", "details": str(e)})

# --- Main execution for direct run or mcp dev ---
if __name__ == "__main__":
    # To run this server with the MCP Inspector:
    # mcp dev chembl_mcp_server.py
    #
    # To install it in Claude Desktop:
    # mcp install chembl_mcp_server.py --name "ChEMBL Tools"

    transport = os.getenv("MCP_TRANSPORT", "stdio")  # Default to stdio, can be overridden by .env or actual env var

    if transport == "stdio":
        print(f"Running ChEMBL MCP Server with stdio transport")
        mcp.run(transport="stdio")
    elif transport == "sse":
        print(f"Running ChEMBL MCP Server with SSE transport on host {mcp.host}, port {mcp.port}")
        mcp.run(transport="sse")
    else:
        raise ValueError(f"Unknown transport: {transport}")

