import sys
import os

# Add the project root to sys.path to allow importing 'engine'
# Assuming this file is in backend/app/routers/
# Adjust the path to go up three levels to the project root (valency/)
# then into the engine directory.
project_root_for_engine = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
# Check if 'engine' directory is directly under project_root_for_engine or if valency is the project_root
# Based on the structure, valency is the root, and engine is a direct child.
# So, we need to add valency (project_root_for_engine) to sys.path
if project_root_for_engine not in sys.path:
    sys.path.insert(0, project_root_for_engine)

from fastapi import APIRouter, HTTPException, Body
# Removed direct import of fastmcp.client.Client
# Assuming 'engine' is installed or accessible in PYTHONPATH
from engine.MCP.MCPAdapter import MCPAdapter 
import os
import json
from typing import Any, List, Optional, Dict
from pydantic import BaseModel, Field

# --- Helper Function for Running Tools ---
async def _run_mcp_tool_handler(
    server_url: str, 
    tool_name: str, 
    tool_args: dict | None = None, 
    request_body: Any = None 
):
    """
    Generic handler to run a tool using MCPAdapter (which internally uses fastmcp.client.Client).
    """
    actual_tool_args = tool_args if tool_args is not None else (request_body if isinstance(request_body, dict) else {})
    
    adapter = MCPAdapter(server_url) # MCPAdapter's __init__ handles /sse in server_url

    try:
        # print(f"Router using MCPAdapter to run tool: {tool_name} with args: {actual_tool_args} on server {server_url}")
        # MCPAdapter.run_tool is expected to return a string (JSON or other format)
        result_content_str = await adapter.run_tool(tool_name, actual_tool_args)
        
        try:
            # Attempt to parse if the string is JSON
            parsed_result = json.loads(result_content_str)
            return {"status": "success", "tool_name": tool_name, "result": parsed_result}
        except json.JSONDecodeError:
            # If not JSON, return the raw string content
            return {"status": "success", "tool_name": tool_name, "result": result_content_str}

    except ConnectionError as ce:
        # print(f"Router error (via MCPAdapter): ConnectionError - {str(ce)}")
        raise HTTPException(status_code=503, detail=f"MCP Adapter connection error for {tool_name} at {server_url}: {str(ce)}")
    except ValueError as ve: 
        # print(f"Router error (via MCPAdapter): ValueError - {str(ve)}")
        raise HTTPException(status_code=400, detail=f"MCP Adapter tool execution error for {tool_name}: {str(ve)}")
    except Exception as e:
        # print(f"Router error (via MCPAdapter): Unexpected Exception - {str(e)}")
        raise HTTPException(status_code=500, detail=f"An unexpected error occurred via MCP Adapter while running {tool_name}: {str(e)}")

# --- ADMET MCP Server Configuration ---
ADMET_MCP_HOST = os.getenv("ADMET_MCP_HOST", "localhost") # Default from ADMET_mcp_server.py was 0.0.0.0, using localhost for client
ADMET_MCP_PORT = os.getenv("ADMET_MCP_PORT", "8055")     # Default from ADMET_mcp_server.py
ADMET_MCP_SERVER_URL = f"http://{ADMET_MCP_HOST}:{ADMET_MCP_PORT}/sse" # MCP server is at /sse
# admet_mcp_adapter = MCPAdapter(server_url=ADMET_MCP_SERVER_URL) # Removed
admet_router = APIRouter(prefix="/mcp/admet", tags=["ADMET MCP Tools"])

# Pydantic Models for ADMET tool arguments
class AdmetGetPredictionArgs(BaseModel):
    smiles: str = Field(..., description="The SMILES string representing the molecule.")

# Explicit ADMET tool endpoint
@admet_router.post(
    "/get_admet_prediction",
    summary="Get ADMET Prediction",
    description="Get ADMET (Absorption, Distribution, Metabolism, Excretion, Toxicity) predictions for a given SMILES string. Returns JSON data."
)
async def admet_get_admet_prediction(tool_args: AdmetGetPredictionArgs = Body(...)):
    return await _run_mcp_tool_handler(ADMET_MCP_SERVER_URL, "get_admet_prediction", tool_args.model_dump(exclude_none=True)) # Pass URL

# --- AlphaFold MCP Server Configuration ---
ALPHAFOLD_MCP_HOST = os.getenv("ALPHAFOLD_MCP_HOST", "localhost") # Default from alphafold_mcp_server.py was 0.0.0.0
ALPHAFOLD_MCP_PORT = os.getenv("ALPHAFOLD_MCP_PORT", "8056")     # Default from alphafold_mcp_server.py was 8052, check your .env
ALPHAFOLD_MCP_SERVER_URL = f"http://{ALPHAFOLD_MCP_HOST}:{ALPHAFOLD_MCP_PORT}/sse" # MCP server is at /sse
# alphafold_mcp_adapter = MCPAdapter(server_url=ALPHAFOLD_MCP_SERVER_URL) # Removed
alphafold_router = APIRouter(prefix="/mcp/alphafold", tags=["AlphaFold MCP Tools"])

# Pydantic Models for AlphaFold tool arguments
class AlphafoldGetPredictionArgs(BaseModel):
    qualifier: str = Field(..., description="UniProt accession (e.g., 'Q5VSL9').")
    sequence_checksum: Optional[str] = Field(None, description="Optional CRC64 checksum of the UniProt sequence.")

class AlphafoldGetUniprotSummaryArgs(BaseModel):
    qualifier: str = Field(..., description="UniProtKB accession number (AC), entry name (ID), or CRC64 checksum of the UniProt sequence (e.g., 'Q5VSL9').")

class AlphafoldGetAnnotationsArgs(BaseModel):
    qualifier: str = Field(..., description="UniProt accession (e.g., 'Q5VSL9').")
    annotation_type: str = Field(..., description="Type of annotation (e.g., 'MUTAGEN' for AlphaMissense).")

@alphafold_router.post(
    "/get_alphafold_prediction",
    summary="Get AlphaFold Prediction Models",
    description="Get all AlphaFold models for a UniProt accession. Returns JSON data from AlphaFold DB."
)
async def alphafold_get_alphafold_prediction(tool_args: AlphafoldGetPredictionArgs = Body(...)):
    return await _run_mcp_tool_handler(ALPHAFOLD_MCP_SERVER_URL, "get_alphafold_prediction", tool_args.model_dump(exclude_none=True)) # Pass URL

@alphafold_router.post(
    "/get_uniprot_summary",
    summary="Get UniProt Summary",
    description="Get summary details for a UniProt residue range. Returns JSON data."
)
async def alphafold_get_uniprot_summary(tool_args: AlphafoldGetUniprotSummaryArgs = Body(...)):
    return await _run_mcp_tool_handler(ALPHAFOLD_MCP_SERVER_URL, "get_uniprot_summary", tool_args.model_dump(exclude_none=True)) # Pass URL

@alphafold_router.post(
    "/get_alphafold_annotations",
    summary="Get AlphaFold Annotations",
    description="Get all annotations for a UniProt residue range. Returns JSON data."
)
async def alphafold_get_alphafold_annotations(tool_args: AlphafoldGetAnnotationsArgs = Body(...)):
    return await _run_mcp_tool_handler(ALPHAFOLD_MCP_SERVER_URL, "get_alphafold_annotations", tool_args.model_dump(exclude_none=True)) # Pass URL

# --- ChEMBL MCP Server Configuration ---
CHEMBL_MCP_HOST = os.getenv("CHEMBL_HOST", "localhost")
CHEMBL_MCP_PORT = os.getenv("CHEMBL_PORT", "8051")
CHEMBL_MCP_SERVER_URL = f"http://{CHEMBL_MCP_HOST}:{CHEMBL_MCP_PORT}/sse"
# chembl_mcp_adapter = MCPAdapter(server_url=CHEMBL_MCP_SERVER_URL) # Removed
chembl_router = APIRouter(prefix="/mcp/chembl", tags=["ChEMBL MCP Tools"])

# Pydantic Models for ChEMBL tool arguments
class ChemblGetMoleculeByIdArgs(BaseModel):
    chembl_id: str = Field(..., description="The ChEMBL ID of the molecule (e.g., 'CHEMBL192').")

class ChemblFindMoleculeByPrefNameArgs(BaseModel):
    pref_name: str = Field(..., description="The preferred name to search for (e.g., 'Aspirin').")

class ChemblFindMoleculeBySynonymArgs(BaseModel):
    synonym: str = Field(..., description="The synonym to search for (case-insensitive exact match).")

class ChemblGetMoleculesByIdsArgs(BaseModel):
    chembl_ids: List[str] = Field(..., description="A list of ChEMBL IDs (e.g., ['CHEMBL25', 'CHEMBL192']).")

class ChemblGetMoleculeImageSvgArgs(BaseModel):
    chembl_id: str = Field(..., description="The ChEMBL ID of the molecule.")

class ChemblFindSimilarMoleculesBySmilesArgs(BaseModel):
    smiles: str = Field(..., description="The SMILES string of the query molecule (e.g., 'CCO').")
    similarity_threshold: int = Field(70, description="Minimum similarity percentage (0-100, default 70).")

class ChemblFindSimilarMoleculesByIdArgs(BaseModel):
    chembl_id: str = Field(..., description="The ChEMBL ID of the query molecule.")
    similarity_threshold: int = Field(70, description="Minimum similarity percentage (0-100, default 70).")

class ChemblGetApprovedDrugsArgs(BaseModel):
    order_by_mw: bool = Field(False, description="If True, sorts results by molecular weight.")

class ChemblGetActivitiesForTargetArgs(BaseModel):
    biological_target_name: str = Field(..., description="Biological target name (e.g., 'hERG').")
    standard_type: Optional[str] = Field("", description="Bioactivity measurement type (e.g., 'IC50', 'Ki'). None for all.")

class ChemblGetActivitiesForMoleculeArgs(BaseModel):
    molecule_chembl_id: str = Field(..., description="ChEMBL ID of the molecule.")
    pchembl_value_exists: bool = Field(True, description="If True, only returns activities with a pChEMBL value.")

class ChemblFindTargetByGeneNameArgs(BaseModel):
    gene_name: str = Field(..., description="Gene name or symbol (e.g., 'BRCA1'). Case-insensitive contains match.")

class ChemblGetMoleculesByFilterArgs(BaseModel):
    filters: Dict[str, Any] = Field(..., description="Django-style filters (e.g., {'molecule_properties__mw_freebase__gte': 200}).")
    order_by: Optional[List[str]] = Field(None, description="Fields to sort by (e.g., ['molecule_properties__mw_freebase']).")

class ChemblGetActivitiesByFilterArgs(BaseModel):
    filters: Dict[str, Any] = Field(..., description="Filters for activity fields (e.g., {'pchembl_value__gte': 6.0}).")
    order_by: Optional[List[str]] = Field(None, description="Fields to sort by (e.g., ['-pchembl_value']).")

class ChemblGetTargetsByFilterArgs(BaseModel):
    filters: Dict[str, Any] = Field(..., description="Filters for target fields (e.g., {'target_type': 'SINGLE PROTEIN'}).")
    order_by: Optional[List[str]] = Field(None, description="Fields to sort by (e.g., ['pref_name']).")

class ChemblSmilesToCtabArgs(BaseModel):
    smiles: str = Field(..., description="SMILES string to convert.")

class ChemblComputeMolecularDescriptorsArgs(BaseModel):
    smiles: str = Field(..., description="SMILES string for descriptor calculation.")

class ChemblComputeStructuralAlertsArgs(BaseModel):
    smiles: str = Field(..., description="SMILES string for structural alert analysis.")

class ChemblStandardizeMoleculeFromSmilesArgs(BaseModel):
    smiles: str = Field(..., description="SMILES string of the molecule to standardize.")

class ChemblGetParentMoleculeFromSmilesArgs(BaseModel):
    smiles: str = Field(..., description="SMILES string of the molecule (possibly salt/mixture).")

# Explicit ChEMBL tool endpoints
@chembl_router.post("/get_molecule_by_chembl_id", summary="Get Molecule by ChEMBL ID", description="Retrieve a specific molecule by its unique ChEMBL ID.")
async def chembl_get_molecule_by_chembl_id(tool_args: ChemblGetMoleculeByIdArgs = Body(...)):
    return await _run_mcp_tool_handler(CHEMBL_MCP_SERVER_URL, "get_molecule_by_chembl_id", tool_args.model_dump(exclude_none=True)) # Pass URL

@chembl_router.post("/find_molecule_by_pref_name", summary="Find Molecule by Preferred Name", description="Search for molecules by their preferred name.")
async def chembl_find_molecule_by_pref_name(tool_args: ChemblFindMoleculeByPrefNameArgs = Body(...)):
    return await _run_mcp_tool_handler(CHEMBL_MCP_SERVER_URL, "find_molecule_by_pref_name", tool_args.model_dump(exclude_none=True)) # Pass URL

@chembl_router.post("/find_molecule_by_synonym", summary="Find Molecule by Synonym", description="Find molecules by their synonyms.")
async def chembl_find_molecule_by_synonym(tool_args: ChemblFindMoleculeBySynonymArgs = Body(...)):
    return await _run_mcp_tool_handler(CHEMBL_MCP_SERVER_URL, "find_molecule_by_synonym", tool_args.model_dump(exclude_none=True)) # Pass URL

@chembl_router.post("/get_molecules_by_chembl_ids", summary="Get Molecules by ChEMBL IDs", description="Retrieve multiple molecules by a list of ChEMBL IDs.")
async def chembl_get_molecules_by_chembl_ids(tool_args: ChemblGetMoleculesByIdsArgs = Body(...)):
    return await _run_mcp_tool_handler(CHEMBL_MCP_SERVER_URL, "get_molecules_by_chembl_ids", tool_args.model_dump(exclude_none=True)) # Pass URL

@chembl_router.post("/get_molecule_image_svg", summary="Get Molecule Image SVG", description="Get the 2D chemical structure image of a molecule as an SVG string.")
async def chembl_get_molecule_image_svg(tool_args: ChemblGetMoleculeImageSvgArgs = Body(...)):
    return await _run_mcp_tool_handler(CHEMBL_MCP_SERVER_URL, "get_molecule_image_svg", tool_args.model_dump(exclude_none=True)) # Pass URL

@chembl_router.post("/find_similar_molecules_by_smiles", summary="Find Similar Molecules by SMILES", description="Find molecules structurally similar to a given SMILES string.")
async def chembl_find_similar_molecules_by_smiles(tool_args: ChemblFindSimilarMoleculesBySmilesArgs = Body(...)):
    return await _run_mcp_tool_handler(CHEMBL_MCP_SERVER_URL, "find_similar_molecules_by_smiles", tool_args.model_dump(exclude_none=True)) # Pass URL

@chembl_router.post("/find_similar_molecules_by_chembl_id", summary="Find Similar Molecules by ChEMBL ID", description="Find molecules structurally similar to a given ChEMBL ID.")
async def chembl_find_similar_molecules_by_chembl_id(tool_args: ChemblFindSimilarMoleculesByIdArgs = Body(...)):
    return await _run_mcp_tool_handler(CHEMBL_MCP_SERVER_URL, "find_similar_molecules_by_chembl_id", tool_args.model_dump(exclude_none=True)) # Pass URL

@chembl_router.post("/get_approved_drugs", summary="Get Approved Drugs", description="Retrieve all drugs that have reached the maximum clinical trial phase.")
async def chembl_get_approved_drugs(tool_args: ChemblGetApprovedDrugsArgs = Body(...)):
    return await _run_mcp_tool_handler(CHEMBL_MCP_SERVER_URL, "get_approved_drugs", tool_args.model_dump(exclude_none=True)) # Pass URL

@chembl_router.post("/get_activities_for_target", summary="Get Activities for Target", description="Fetch bioactivity data for a specific biological target.")
async def chembl_get_activities_for_target(tool_args: ChemblGetActivitiesForTargetArgs = Body(...)):
    return await _run_mcp_tool_handler(CHEMBL_MCP_SERVER_URL, "get_activities_for_target", tool_args.model_dump(exclude_none=True)) # Pass URL

@chembl_router.post("/get_activities_for_molecule", summary="Get Activities for Molecule", description="Retrieve all recorded bioactivities for a specific molecule.")
async def chembl_get_activities_for_molecule(tool_args: ChemblGetActivitiesForMoleculeArgs = Body(...)):
    return await _run_mcp_tool_handler(CHEMBL_MCP_SERVER_URL, "get_activities_for_molecule", tool_args.model_dump(exclude_none=True)) # Pass URL

@chembl_router.post("/find_target_by_gene_name", summary="Find Target by Gene Name", description="Search for biological targets by a gene name or symbol.")
async def chembl_find_target_by_gene_name(tool_args: ChemblFindTargetByGeneNameArgs = Body(...)):
    return await _run_mcp_tool_handler(CHEMBL_MCP_SERVER_URL, "find_target_by_gene_name", tool_args.model_dump(exclude_none=True)) # Pass URL

@chembl_router.post("/get_molecules_by_filter", summary="Get Molecules by Filter", description="Retrieve molecules based on a custom set of filter conditions.")
async def chembl_get_molecules_by_filter(tool_args: ChemblGetMoleculesByFilterArgs = Body(...)):
    return await _run_mcp_tool_handler(CHEMBL_MCP_SERVER_URL, "get_molecules_by_filter", tool_args.model_dump(exclude_none=True)) # Pass URL

@chembl_router.post("/get_activities_by_filter", summary="Get Activities by Filter", description="Retrieve bioactivity data based on a custom set of filter conditions.")
async def chembl_get_activities_by_filter(tool_args: ChemblGetActivitiesByFilterArgs = Body(...)):
    return await _run_mcp_tool_handler(CHEMBL_MCP_SERVER_URL, "get_activities_by_filter", tool_args.model_dump(exclude_none=True)) # Pass URL

@chembl_router.post("/get_targets_by_filter", summary="Get Targets by Filter", description="Retrieve biological targets based on a custom set of filter conditions.")
async def chembl_get_targets_by_filter(tool_args: ChemblGetTargetsByFilterArgs = Body(...)):
    return await _run_mcp_tool_handler(CHEMBL_MCP_SERVER_URL, "get_targets_by_filter", tool_args.model_dump(exclude_none=True)) # Pass URL

@chembl_router.post("/smiles_to_ctab", summary="SMILES to CTAB", description="Convert a SMILES string to a CTAB block.")
async def chembl_smiles_to_ctab(tool_args: ChemblSmilesToCtabArgs = Body(...)):
    return await _run_mcp_tool_handler(CHEMBL_MCP_SERVER_URL, "smiles_to_ctab", tool_args.model_dump(exclude_none=True)) # Pass URL

@chembl_router.post("/compute_molecular_descriptors", summary="Compute Molecular Descriptors", description="Calculate physicochemical properties for a molecule from its SMILES string.")
async def chembl_compute_molecular_descriptors(tool_args: ChemblComputeMolecularDescriptorsArgs = Body(...)):
    return await _run_mcp_tool_handler(CHEMBL_MCP_SERVER_URL, "compute_molecular_descriptors", tool_args.model_dump(exclude_none=True)) # Pass URL

@chembl_router.post("/compute_structural_alerts", summary="Compute Structural Alerts", description="Identify known structural alerts within a molecule from its SMILES string.")
async def chembl_compute_structural_alerts(tool_args: ChemblComputeStructuralAlertsArgs = Body(...)):
    return await _run_mcp_tool_handler(CHEMBL_MCP_SERVER_URL, "compute_structural_alerts", tool_args.model_dump(exclude_none=True)) # Pass URL

@chembl_router.post("/standardize_molecule_from_smiles", summary="Standardize Molecule from SMILES", description="Standardize a molecular structure provided as a SMILES string.")
async def chembl_standardize_molecule_from_smiles(tool_args: ChemblStandardizeMoleculeFromSmilesArgs = Body(...)):
    return await _run_mcp_tool_handler(CHEMBL_MCP_SERVER_URL, "standardize_molecule_from_smiles", tool_args.model_dump(exclude_none=True)) # Pass URL

@chembl_router.post("/get_parent_molecule_from_smiles", summary="Get Parent Molecule from SMILES", description="Identify and return the parent structure for a molecule (e.g., removing counter-ions).")
async def chembl_get_parent_molecule_from_smiles(tool_args: ChemblGetParentMoleculeFromSmilesArgs = Body(...)):
    return await _run_mcp_tool_handler(CHEMBL_MCP_SERVER_URL, "get_parent_molecule_from_smiles", tool_args.model_dump(exclude_none=True)) # Pass URL


# --- PubChem MCP Server Configuration ---
PUBCHEM_MCP_HOST = os.getenv("PUBCHEM_MCP_HOST", "localhost")
PUBCHEM_MCP_PORT = os.getenv("PUBCHEM_MCP_PORT", "8054")
PUBCHEM_MCP_SERVER_URL = f"http://{PUBCHEM_MCP_HOST}:{PUBCHEM_MCP_PORT}/sse"
# pubchem_mcp_adapter = MCPAdapter(server_url=PUBCHEM_MCP_SERVER_URL) # Removed
pubchem_router = APIRouter(prefix="/mcp/pubchem", tags=["PubChem MCP Tools"])

# Pydantic Models for PubChem tool arguments
# Note: These models are based on inferred arguments due to limited context on exact tool signatures from pubchem_mcp_server.py.
class PubchemGetCompoundByCidArgs(BaseModel):
    cid: str = Field(..., description="PubChem Compound ID (CID).")

class PubchemGetCidsByNameArgs(BaseModel):
    name: str = Field(..., description="Name of the compound to search for.")
    match_type: Optional[str] = Field(None, description="Type of name match (e.g., 'exact', 'contains'). Refer to PubChem PUG REST for details.")

class PubchemGetCompoundPropertiesArgs(BaseModel):
    cid: str = Field(..., description="PubChem Compound ID (CID).")
    properties: Optional[List[str]] = Field(None, description="Comma-separated string or list of properties to retrieve. Refer to PubChem PUG REST for available properties.")

class PubchemGetCompoundSynonymsByCidArgs(BaseModel):
    cid: str = Field(..., description="PubChem Compound ID (CID).")

class PubchemGetCompoundImagePubchemUrlArgs(BaseModel):
    cid: str = Field(..., description="PubChem Compound ID (CID).")
    image_size: Optional[str] = Field("small", description="Desired image size (e.g., 'small', 'large').") # Assuming, check PubChem docs

class PubchemGetCidsBySmilesArgs(BaseModel):
    smiles: str = Field(..., description="SMILES string.")
    search_type: Optional[str] = Field(None, description="Type of SMILES search (e.g., 'similarity', 'substructure').")

class PubchemGetCidsByInchikeyArgs(BaseModel):
    inchikey: str = Field(..., description="InChIKey.")

class PubchemFastIdentitySearchByCidArgs(BaseModel):
    cid: str = Field(..., description="PubChem Compound ID (CID) for query.")
    threshold: Optional[int] = Field(None, description="Identity threshold for the search.") # e.g. 90 for 90%

class PubchemFastSubstructureSearchBySmilesArgs(BaseModel):
    smiles: str = Field(..., description="SMILES string for substructure query.")

class PubchemFastSimilarity2dSearchByCidArgs(BaseModel):
    cid: str = Field(..., description="PubChem Compound ID (CID) for similarity query.")
    threshold: Optional[int] = Field(None, description="Tanimoto similarity threshold (e.g., 80 for 80%).")

class PubchemGetCidsByXrefArgs(BaseModel):
    xref: str = Field(..., description="Cross-reference ID.")
    namespace: str = Field(..., description="Namespace of the cross-reference (e.g., 'RegistryID', 'RN', 'PubMedID').")

class PubchemGetCidsByMassArgs(BaseModel):
    mass: float = Field(..., description="Molecular mass for search.")
    tolerance: Optional[float] = Field(None, description="Mass tolerance for the search.")

# Explicit PubChem tool endpoints
@pubchem_router.post("/get_compound_by_cid", summary="Get Compound by CID", description="Retrieve compound details from PubChem using its Compound ID (CID).")
async def pubchem_get_compound_by_cid(tool_args: PubchemGetCompoundByCidArgs = Body(...)):
    return await _run_mcp_tool_handler(PUBCHEM_MCP_SERVER_URL, "get_compound_by_cid", tool_args.model_dump(exclude_none=True)) # Pass URL

@pubchem_router.post("/get_cids_by_name", summary="Get CIDs by Name", description="Search for PubChem Compound IDs (CIDs) by compound name.")
async def pubchem_get_cids_by_name(tool_args: PubchemGetCidsByNameArgs = Body(...)):
    return await _run_mcp_tool_handler(PUBCHEM_MCP_SERVER_URL, "get_cids_by_name", tool_args.model_dump(exclude_none=True)) # Pass URL

@pubchem_router.post("/get_compound_properties", summary="Get Compound Properties", description="Retrieve specific properties for a given PubChem CID.")
async def pubchem_get_compound_properties(tool_args: PubchemGetCompoundPropertiesArgs = Body(...)):
    return await _run_mcp_tool_handler(PUBCHEM_MCP_SERVER_URL, "get_compound_properties", tool_args.model_dump(exclude_none=True)) # Pass URL

@pubchem_router.post("/get_compound_synonyms_by_cid", summary="Get Compound Synonyms by CID", description="Retrieve synonyms for a compound using its PubChem CID.")
async def pubchem_get_compound_synonyms_by_cid(tool_args: PubchemGetCompoundSynonymsByCidArgs = Body(...)):
    return await _run_mcp_tool_handler(PUBCHEM_MCP_SERVER_URL, "get_compound_synonyms_by_cid", tool_args.model_dump(exclude_none=True)) # Pass URL

@pubchem_router.post("/get_compound_image_pubchem_url", summary="Get Compound Image URL", description="Get a URL for the image of a compound from PubChem using its CID.")
async def pubchem_get_compound_image_pubchem_url(tool_args: PubchemGetCompoundImagePubchemUrlArgs = Body(...)):
    return await _run_mcp_tool_handler(PUBCHEM_MCP_SERVER_URL, "get_compound_image_pubchem_url", tool_args.model_dump(exclude_none=True)) # Pass URL

@pubchem_router.post("/get_cids_by_smiles", summary="Get CIDs by SMILES", description="Search for PubChem CIDs using a SMILES string.")
async def pubchem_get_cids_by_smiles(tool_args: PubchemGetCidsBySmilesArgs = Body(...)):
    return await _run_mcp_tool_handler(PUBCHEM_MCP_SERVER_URL, "get_cids_by_smiles", tool_args.model_dump(exclude_none=True)) # Pass URL

@pubchem_router.post("/get_cids_by_inchikey", summary="Get CIDs by InChIKey", description="Search for PubChem CIDs using an InChIKey.")
async def pubchem_get_cids_by_inchikey(tool_args: PubchemGetCidsByInchikeyArgs = Body(...)):
    return await _run_mcp_tool_handler(PUBCHEM_MCP_SERVER_URL, "get_cids_by_inchikey", tool_args.model_dump(exclude_none=True)) # Pass URL

@pubchem_router.post("/fast_identity_search_by_cid", summary="Fast Identity Search by CID", description="Perform a fast identity search for compounds similar to a given PubChem CID.")
async def pubchem_fast_identity_search_by_cid(tool_args: PubchemFastIdentitySearchByCidArgs = Body(...)):
    return await _run_mcp_tool_handler(PUBCHEM_MCP_SERVER_URL, "fast_identity_search_by_cid", tool_args.model_dump(exclude_none=True)) # Pass URL

@pubchem_router.post("/fast_substructure_search_by_smiles", summary="Fast Substructure Search by SMILES", description="Perform a fast substructure search using a SMILES string.")
async def pubchem_fast_substructure_search_by_smiles(tool_args: PubchemFastSubstructureSearchBySmilesArgs = Body(...)):
    return await _run_mcp_tool_handler(PUBCHEM_MCP_SERVER_URL, "fast_substructure_search_by_smiles", tool_args.model_dump(exclude_none=True)) # Pass URL

@pubchem_router.post("/fast_similarity_2d_search_by_cid", summary="Fast 2D Similarity Search by CID", description="Perform a fast 2D similarity search for compounds similar to a given PubChem CID.")
async def pubchem_fast_similarity_2d_search_by_cid(tool_args: PubchemFastSimilarity2dSearchByCidArgs = Body(...)):
    return await _run_mcp_tool_handler(PUBCHEM_MCP_SERVER_URL, "fast_similarity_2d_search_by_cid", tool_args.model_dump(exclude_none=True)) # Pass URL

@pubchem_router.post("/get_cids_by_xref", summary="Get CIDs by Cross-Reference (XRef)", description="Search for PubChem CIDs using a cross-reference ID and its namespace.")
async def pubchem_get_cids_by_xref(tool_args: PubchemGetCidsByXrefArgs = Body(...)):
    return await _run_mcp_tool_handler(PUBCHEM_MCP_SERVER_URL, "get_cids_by_xref", tool_args.model_dump(exclude_none=True)) # Pass URL

@pubchem_router.post("/get_cids_by_mass", summary="Get CIDs by Mass", description="Search for PubChem CIDs using molecular mass and an optional tolerance.")
async def pubchem_get_cids_by_mass(tool_args: PubchemGetCidsByMassArgs = Body(...)):
    return await _run_mcp_tool_handler(PUBCHEM_MCP_SERVER_URL, "get_cids_by_mass", tool_args.model_dump(exclude_none=True)) # Pass URL

# --- RXN for Chemistry MCP Server Configuration ---
RXN_MCP_HOST = os.getenv("RXN_MCP_HOST", "localhost") # Default from rxn_mcp_server.py was 0.0.0.0
RXN_MCP_PORT = os.getenv("RXN_MCP_PORT", "8057")     # Default from rxn_mcp_server.py
RXN_MCP_SERVER_URL = f"http://{RXN_MCP_HOST}:{RXN_MCP_PORT}/sse"
rxn_router = APIRouter(prefix="/mcp/rxn", tags=["RXN for Chemistry MCP Tools"])

# Pydantic Models for RXN tool arguments
class RxnPredictReactionSmilesArgs(BaseModel):
    precursors_smiles: str = Field(..., description="SMILES of precursor molecules, separated by a dot (e.g., 'CCO.Cl').")

class RxnGetReactionPredictionResultsArgs(BaseModel):
    prediction_id: str = Field(..., description="ID of the prediction task.")

class RxnPredictAutomaticRetrosynthesisSmilesArgs(BaseModel):
    product_smiles: str = Field(..., description="SMILES string of the target molecule.")
    availability_pricing_threshold: float = Field(0.0, description="Exclude precursors with price higher than this (0.0 means no filtering).")
    available_smiles_list: Optional[List[str]] = Field(None, description="List of SMILES that are considered available.")
    exclude_smiles_list: Optional[List[str]] = Field(None, description="List of SMILES to exclude from precursors.")
    exclude_substructures_smiles_list: Optional[List[str]] = Field(None, description="List of SMILES substructures to exclude from precursors.")
    fap_value: float = Field(0.6, description="Forward Acceptance Probability threshold (0.0 to 1.0).")
    max_steps: int = Field(3, description="Maximum number of retrosynthetic steps.")
    n_beams: int = Field(3, description="Number of beams to use in the search.")
    reaction_class_priority_list: Optional[List[str]] = Field(None, description="List of reaction classes to prioritize.")
    project_id: Optional[str] = Field(None, description="Optional ID of the project to associate this prediction with.")

class RxnGetRetrosynthesisPredictionResultsArgs(BaseModel):
    prediction_id: str = Field(..., description="ID of the retrosynthesis prediction task.")
    project_id: Optional[str] = Field(None, description="Optional ID of the project this prediction belongs to.")

class RxnPredictReagentsSmilesArgs(BaseModel):
    starting_material_smiles: str = Field(..., description="SMILES string of the starting material.")
    product_smiles: str = Field(..., description="SMILES string of the product.")
    project_id: Optional[str] = Field(None, description="Optional ID of the project to associate this prediction with.")

class RxnGetReagentPredictionResultsArgs(BaseModel):
    prediction_id: str = Field(..., description="ID of the reagent prediction task.")
    project_id: Optional[str] = Field(None, description="Optional ID of the project this prediction belongs to.")

class RxnGetAtomMappingForReactionSmilesArgs(BaseModel):
    reaction_smiles_list: List[str] = Field(..., description="List of reaction SMILES strings (e.g., ['CCO.Cl>>CCCl.O']).")

class RxnTranslateParagraphToActionsArgs(BaseModel):
    paragraph_text: str = Field(..., description="Text paragraph describing the chemical procedure.")

class RxnDigitizeReactionFromFileIdArgs(BaseModel):
    file_id: str = Field(..., description="ID of the file uploaded to RXN.")
    project_id: Optional[str] = Field(None, description="Optional ID of the project to associate this digitization with.")

class RxnUploadFileArgs(BaseModel):
    server_file_path: str = Field(..., description="Absolute path to the file on the server where the MCP is running.")

class RxnPredictReactionBatchSmilesArgs(BaseModel):
    precursors_smiles_list: List[str] = Field(..., description="List of precursor SMILES strings for batch prediction.")

class RxnGetReactionBatchPredictionResultsArgs(BaseModel):
    task_id: str = Field(..., description="ID of the batch prediction task.")

class RxnPredictReactionBatchTopnSmilesArgs(BaseModel):
    precursors_smiles_lists: List[List[str]] = Field(..., description="List of lists of precursor SMILES for top-N batch prediction.")
    top_n: int = Field(..., description="Number of top predictions to return for each reaction.")

class RxnGetReactionBatchTopnPredictionResultsArgs(BaseModel):
    task_id: str = Field(..., description="ID of the batch top-N prediction task.")

class RxnCreateSynthesisFromSequenceArgs(BaseModel):
    sequence_id: str = Field(..., description="ID of the retrosynthetic sequence.")
    project_id: Optional[str] = Field(None, description="Optional ID of the project.")

class RxnGetSynthesisNodeIdsArgs(BaseModel):
    synthesis_id: str = Field(..., description="ID of the synthesis plan.")
    project_id: Optional[str] = Field(None, description="Optional ID of the project.")

class RxnGetSynthesisNodeReactionSettingsArgs(BaseModel):
    synthesis_id: str = Field(..., description="ID of the synthesis plan.")
    node_id: str = Field(..., description="ID of the node within the synthesis plan.")
    project_id: Optional[str] = Field(None, description="Optional ID of the project.")

class RxnUpdateSynthesisNodeReactionSettingsArgs(BaseModel):
    synthesis_id: str = Field(..., description="ID of the synthesis plan.")
    node_id: str = Field(..., description="ID of the node within the synthesis plan.")
    actions: List[Dict[str, Any]] = Field(..., description="List of action dictionaries defining the procedure.")
    product_smiles: str = Field(..., description="SMILES string of the product for this reaction step.")
    project_id: Optional[str] = Field(None, description="Optional ID of the project.")

class RxnCreateProjectArgs(BaseModel):
    name: str = Field(..., description="Name for the new project.")

class RxnSetCurrentProjectArgs(BaseModel):
    project_id: str = Field(..., description="ID of the project to set as active.")

# No args for rxn_get_current_project_id

# Explicit RXN tool endpoints
@rxn_router.post("/predict_reaction_smiles", summary="Predict Reaction from SMILES", description="Predicts the product of a chemical reaction given precursors as SMILES.")
async def rxn_predict_reaction_smiles(tool_args: RxnPredictReactionSmilesArgs = Body(...)):
    return await _run_mcp_tool_handler(RXN_MCP_SERVER_URL, "predict_reaction_smiles", tool_args.model_dump(exclude_none=True))

@rxn_router.post("/get_reaction_prediction_results", summary="Get Reaction Prediction Results", description="Retrieves the results of a reaction prediction task.")
async def rxn_get_reaction_prediction_results(tool_args: RxnGetReactionPredictionResultsArgs = Body(...)):
    return await _run_mcp_tool_handler(RXN_MCP_SERVER_URL, "get_reaction_prediction_results", tool_args.model_dump(exclude_none=True))

@rxn_router.post("/predict_automatic_retrosynthesis_smiles", summary="Predict Automatic Retrosynthesis from SMILES", description="Predicts possible retrosynthetic routes for a target molecule SMILES.")
async def rxn_predict_automatic_retrosynthesis_smiles(tool_args: RxnPredictAutomaticRetrosynthesisSmilesArgs = Body(...)):
    return await _run_mcp_tool_handler(RXN_MCP_SERVER_URL, "predict_automatic_retrosynthesis_smiles", tool_args.model_dump(exclude_none=True))

@rxn_router.post("/get_retrosynthesis_prediction_results", summary="Get Retrosynthesis Prediction Results", description="Retrieves the results of an automatic retrosynthesis prediction task.")
async def rxn_get_retrosynthesis_prediction_results(tool_args: RxnGetRetrosynthesisPredictionResultsArgs = Body(...)):
    return await _run_mcp_tool_handler(RXN_MCP_SERVER_URL, "get_retrosynthesis_prediction_results", tool_args.model_dump(exclude_none=True))

@rxn_router.post("/predict_reagents_smiles", summary="Predict Reagents from SMILES", description="Predicts reagents needed to convert a starting material to a product.")
async def rxn_predict_reagents_smiles(tool_args: RxnPredictReagentsSmilesArgs = Body(...)):
    return await _run_mcp_tool_handler(RXN_MCP_SERVER_URL, "predict_reagents_smiles", tool_args.model_dump(exclude_none=True))

@rxn_router.post("/get_reagent_prediction_results", summary="Get Reagent Prediction Results", description="Retrieves the results of a reagent prediction task.")
async def rxn_get_reagent_prediction_results(tool_args: RxnGetReagentPredictionResultsArgs = Body(...)):
    return await _run_mcp_tool_handler(RXN_MCP_SERVER_URL, "get_reagent_prediction_results", tool_args.model_dump(exclude_none=True))

@rxn_router.post("/get_atom_mapping_for_reaction_smiles", summary="Get Atom Mapping for Reaction SMILES", description="Performs atom mapping for a list of chemical reactions.")
async def rxn_get_atom_mapping_for_reaction_smiles(tool_args: RxnGetAtomMappingForReactionSmilesArgs = Body(...)):
    return await _run_mcp_tool_handler(RXN_MCP_SERVER_URL, "get_atom_mapping_for_reaction_smiles", tool_args.model_dump(exclude_none=True))

@rxn_router.post("/translate_paragraph_to_actions", summary="Translate Paragraph to Actions", description="Translates a natural language chemical procedure into machine-readable actions.")
async def rxn_translate_paragraph_to_actions(tool_args: RxnTranslateParagraphToActionsArgs = Body(...)):
    return await _run_mcp_tool_handler(RXN_MCP_SERVER_URL, "translate_paragraph_to_actions", tool_args.model_dump(exclude_none=True))

@rxn_router.post("/digitize_reaction_from_file_id", summary="Digitize Reaction from File ID", description="Starts digitization for a reaction scheme from an uploaded file ID.")
async def rxn_digitize_reaction_from_file_id(tool_args: RxnDigitizeReactionFromFileIdArgs = Body(...)):
    return await _run_mcp_tool_handler(RXN_MCP_SERVER_URL, "digitize_reaction_from_file_id", tool_args.model_dump(exclude_none=True))

@rxn_router.post("/rxn_upload_file", summary="Upload File to RXN", description="Uploads a file to the RXN server.")
async def rxn_upload_file(tool_args: RxnUploadFileArgs = Body(...)):
    return await _run_mcp_tool_handler(RXN_MCP_SERVER_URL, "rxn_upload_file", tool_args.model_dump(exclude_none=True))

@rxn_router.post("/predict_reaction_batch_smiles", summary="Predict Reaction Batch from SMILES", description="Predicts products for a batch of chemical reactions.")
async def rxn_predict_reaction_batch_smiles(tool_args: RxnPredictReactionBatchSmilesArgs = Body(...)):
    return await _run_mcp_tool_handler(RXN_MCP_SERVER_URL, "predict_reaction_batch_smiles", tool_args.model_dump(exclude_none=True))

@rxn_router.post("/get_reaction_batch_prediction_results", summary="Get Reaction Batch Prediction Results", description="Retrieves results of a batch reaction prediction task.")
async def rxn_get_reaction_batch_prediction_results(tool_args: RxnGetReactionBatchPredictionResultsArgs = Body(...)):
    return await _run_mcp_tool_handler(RXN_MCP_SERVER_URL, "get_reaction_batch_prediction_results", tool_args.model_dump(exclude_none=True))

@rxn_router.post("/predict_reaction_batch_topn_smiles", summary="Predict Reaction Batch Top-N SMILES", description="Predicts top N products for a batch of reactions.")
async def rxn_predict_reaction_batch_topn_smiles(tool_args: RxnPredictReactionBatchTopnSmilesArgs = Body(...)):
    return await _run_mcp_tool_handler(RXN_MCP_SERVER_URL, "predict_reaction_batch_topn_smiles", tool_args.model_dump(exclude_none=True))

@rxn_router.post("/get_reaction_batch_topn_prediction_results", summary="Get Reaction Batch Top-N Prediction Results", description="Retrieves results of a batch top-N reaction prediction task.")
async def rxn_get_reaction_batch_topn_prediction_results(tool_args: RxnGetReactionBatchTopnPredictionResultsArgs = Body(...)):
    return await _run_mcp_tool_handler(RXN_MCP_SERVER_URL, "get_reaction_batch_topn_prediction_results", tool_args.model_dump(exclude_none=True))

@rxn_router.post("/create_synthesis_from_sequence", summary="Create Synthesis from Sequence", description="Creates an RXN synthesis plan from a retrosynthetic sequence ID.")
async def rxn_create_synthesis_from_sequence(tool_args: RxnCreateSynthesisFromSequenceArgs = Body(...)):
    return await _run_mcp_tool_handler(RXN_MCP_SERVER_URL, "create_synthesis_from_sequence", tool_args.model_dump(exclude_none=True))

@rxn_router.post("/get_synthesis_node_ids", summary="Get Synthesis Node IDs", description="Retrieves all node IDs for a given synthesis plan.")
async def rxn_get_synthesis_node_ids(tool_args: RxnGetSynthesisNodeIdsArgs = Body(...)):
    return await _run_mcp_tool_handler(RXN_MCP_SERVER_URL, "get_synthesis_node_ids", tool_args.model_dump(exclude_none=True))

@rxn_router.post("/get_synthesis_node_reaction_settings", summary="Get Synthesis Node Reaction Settings", description="Retrieves reaction settings for a specific node in a synthesis plan.")
async def rxn_get_synthesis_node_reaction_settings(tool_args: RxnGetSynthesisNodeReactionSettingsArgs = Body(...)):
    return await _run_mcp_tool_handler(RXN_MCP_SERVER_URL, "get_synthesis_node_reaction_settings", tool_args.model_dump(exclude_none=True))

@rxn_router.post("/update_synthesis_node_reaction_settings", summary="Update Synthesis Node Reaction Settings", description="Updates reaction settings for a specific node in a synthesis plan.")
async def rxn_update_synthesis_node_reaction_settings(tool_args: RxnUpdateSynthesisNodeReactionSettingsArgs = Body(...)):
    return await _run_mcp_tool_handler(RXN_MCP_SERVER_URL, "update_synthesis_node_reaction_settings", tool_args.model_dump(exclude_none=True))

@rxn_router.post("/rxn_create_project", summary="Create RXN Project", description="Creates a new project in RXN for Chemistry.")
async def rxn_create_project(tool_args: RxnCreateProjectArgs = Body(...)):
    return await _run_mcp_tool_handler(RXN_MCP_SERVER_URL, "rxn_create_project", tool_args.model_dump(exclude_none=True))

@rxn_router.post("/rxn_set_current_project", summary="Set Current RXN Project", description="Sets the current active project for the RXN wrapper instance.")
async def rxn_set_current_project(tool_args: RxnSetCurrentProjectArgs = Body(...)):
    return await _run_mcp_tool_handler(RXN_MCP_SERVER_URL, "rxn_set_current_project", tool_args.model_dump(exclude_none=True))

@rxn_router.post("/rxn_get_current_project_id", summary="Get Current RXN Project ID", description="Gets the ID of the currently active project in the RXN wrapper.")
async def rxn_get_current_project_id(): # No arguments
    return await _run_mcp_tool_handler(RXN_MCP_SERVER_URL, "rxn_get_current_project_id", {})

# --- BRICS MCP Server Configuration ---
BRICS_MCP_HOST = os.getenv("BRICS_MCP_HOST", "localhost") # Default from brics_mcp_server.py was 0.0.0.0
BRICS_MCP_PORT = os.getenv("BRICS_MCP_PORT", "8058")     # Default from brics_mcp_server.py
BRICS_MCP_SERVER_URL = f"http://{BRICS_MCP_HOST}:{BRICS_MCP_PORT}/sse"
brics_router = APIRouter(prefix="/mcp/brics", tags=["BRICS MCP Tools"])

# Pydantic Models for BRICS tool arguments
class BricsGetCandidatesArgs(BaseModel):
    smiles_list: List[str] = Field(..., description="A list of SMILES strings.")
    is_polymer: bool = Field(False, description="A boolean flag indicating if the input SMILES are polymers.")

# Explicit BRICS tool endpoint
@brics_router.post(
    "/get_brics_candidates",
    summary="Get BRICS Candidates",
    description="Generate molecular candidates from a list of SMILES strings using the BRICS algorithm."
)
async def brics_get_brics_candidates(tool_args: BricsGetCandidatesArgs = Body(...)):
    return await _run_mcp_tool_handler(BRICS_MCP_SERVER_URL, "get_brics_candidates", tool_args.model_dump(exclude_none=True))

# --- List of all routers to be included in the main FastAPI app ---
# In your main.py or equivalent FastAPI app setup:
# from .mcp_adapter_router import all_mcp_routers
# for router in all_mcp_routers:
#     app.include_router(router)


# The original generic router is now superseded by specific ones.
# If you still need a generic one, it could be adapted or kept.
# For example, the old router code:
# router = APIRouter()
# DEFAULT_MCP_SERVER_URL = "http://localhost:8055/sse"
# MCP_SERVER_HOST = os.getenv("CHEMICAL_EXPLORER_MCP_HOST", "localhost")
# MCP_SERVER_PORT = os.getenv("CHEMICAL_EXPLORER_MCP_PORT", "8055")
# MCP_SERVER_URL = f"http://{MCP_SERVER_HOST}:{MCP_SERVER_PORT}/sse"
# mcp_adapter = MCPAdapter(server_url=MCP_SERVER_URL)
# @router.post("/run_tool/")
# async def run_mcp_tool(tool_name: str = Body(...), tool_args: dict = Body(...)):
#     return await _run_mcp_tool_handler(mcp_adapter, tool_name, tool_args)
# all_mcp_routers.append(router) # If you want to keep it

# Pydantic Models for RCSB PDB tool arguments
class RcsbPdbTextSearchArgs(BaseModel):
    query_string: str = Field(..., description="The text to search for (e.g., \"hemoglobin\").")
    return_type: str = Field("entry", description="Type of identifiers to return (e.g., 'entry', 'polymer_entity').")
    results_verbosity: str = Field("compact", description="Verbosity of results ('compact', 'minimal', 'verbose').")

class RcsbPdbAttributeSearchArgs(BaseModel):
    attribute_path: str = Field(..., description="Path of the attribute to query (e.g., 'exptl.method').")
    operator: str = Field(..., description="Comparison operator (e.g., 'exact_match', 'greater').")
    value: Any = Field(..., description="Value to compare against. For 'in' operator, this should be a list.")
    return_type: str = Field("entry", description="Type of identifiers to return.")
    results_verbosity: str = Field("compact", description="Verbosity of results.")

class RcsbPdbCombinedSearchArgs(BaseModel):
    text_query_string: str = Field(..., description="Main text query.")
    attribute_filters: List[Dict[str, Any]] = Field(..., description='''List of attribute filter dicts. Each dict: {"attribute_path": "...", "operator": "...", "value": ...}. Example: [{"attribute_path": "rcsb_struct_symmetry.symbol", "operator": "exact_match", "value": "C2"}]''')
    logical_operator: str = Field("and", description="How to combine text and attribute filters ('and' or 'or').")
    return_type: str = Field("entry", description="Type of identifiers to return.")
    results_verbosity: str = Field("compact", description="Verbosity of results.")

class RcsbPdbSequenceIdentitySearchArgs(BaseModel):
    sequence: str = Field(..., description="Protein, DNA, or RNA sequence string.")
    identity_cutoff: float = Field(0.9, description="Minimum sequence identity (0.0 to 1.0, default 0.9).")
    e_value_cutoff: float = Field(1.0, description="Maximum E-value for the match (default 1.0).")
    sequence_type: str = Field("protein", description="Type of sequence ('protein', 'dna', 'rna', default 'protein').")
    return_type: str = Field("polymer_entity", description="Type of identifiers to return (default 'polymer_entity').")

class RcsbPdbSequenceMotifSearchArgs(BaseModel):
    motif_pattern: str = Field(..., description='''Motif pattern (e.g., "C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H." for PROSITE).''')
    pattern_type: str = Field("prosite", description="Type of pattern ('prosite', 'regex', 'simple', default 'prosite').")
    sequence_type: str = Field("protein", description="Type of sequence (default 'protein').")
    return_type: str = Field("polymer_entity", description="Type of identifiers to return (default 'polymer_entity').")

class RcsbPdbStructSimilarityEntryIdArgs(BaseModel):
    entry_id: str = Field(..., description="PDB ID of the query structure (e.g., '4HHB').")
    assembly_id: str = Field("1", description="Assembly ID of the query structure (default '1').")
    operator: str = Field("strict_shape_match", description="Similarity operator ('strict_shape_match' or 'relaxed_shape_match').")
    target_search_space: str = Field("assembly", description="What to compare against ('assembly' or 'polymer_entity_instance').")
    return_type: str = Field("assembly", description="Type of identifiers to return.")

class RcsbPdbStructSimilarityFileUrlArgs(BaseModel):
    file_url: str = Field(..., description="URL to the structure file (e.g., 'https://files.rcsb.org/view/4HHB.cif').")
    file_format: str = Field(..., description="Format of the file ('cif', 'bcif', 'pdb', 'cif.gz', 'pdb.gz').")
    operator: str = Field("strict_shape_match", description="Similarity operator.")
    target_search_space: str = Field("assembly", description="What to compare against.")
    return_type: str = Field("assembly", description="Type of identifiers to return.")

class RcsbPdbStructMotifEntryIdArgs(BaseModel):
    entry_id: str = Field(..., description="PDB ID defining the motif (e.g., '2MNR').")
    residues: List[Dict[str, Any]] = Field(..., description='''List of residue definitions. Each dict: {"chain_id": "A", "struct_oper_id": "1", "label_seq_id": 192, "exchanges": ["LYS", "HIS"]}''')
    backbone_distance_tolerance: int = Field(1, description="Allowed backbone distance tolerance in Angstrom (0-3, default 1).")
    side_chain_distance_tolerance: int = Field(1, description="Allowed side-chain distance tolerance in Angstrom (0-3, default 1).")
    angle_tolerance: int = Field(1, description="Allowed angle tolerance in multiples of 20 degrees (0-3, default 1).")
    rmsd_cutoff: float = Field(2.0, description="Threshold above which hits will be filtered by RMSD (>=0, default 2.0).")
    return_type: str = Field("polymer_entity", description="Type of identifiers to return (default 'polymer_entity').")

class RcsbPdbGetTermFacetsArgs(BaseModel):
    attribute_query_dict: Dict[str, Any] = Field(..., description='''Base attribute query dict (e.g., {"attribute_path": "rcsb_entry_info.structure_determination_methodology", "operator": "exact_match", "value": "experimental"}).''')
    facet_name: str = Field(..., description="Descriptive name for this facet (e.g., 'Experimental Methods').")
    facet_attribute: str = Field(..., description="Attribute to aggregate on (e.g., 'exptl.method').")
    min_interval_population: int = Field(1, description="Minimum count for a bucket to be returned (default 1).")
    max_num_intervals: int = Field(10, description="Maximum number of buckets to return (default 10).")

class RcsbPdbGetHistogramFacetsArgs(BaseModel):
    attribute_query_dict: Dict[str, Any] = Field(..., description="Dictionary for the base attribute query.")
    facet_name: str = Field(..., description="Descriptive name for the facet (e.g., 'Formula Weight Distribution').")
    facet_attribute: str = Field(..., description="Numeric attribute for histogram (e.g., 'rcsb_polymer_entity.formula_weight').")
    interval: float = Field(..., description="Size of the histogram intervals.")
    min_interval_population: int = Field(1, description="Minimum count for a bucket (default 1).")
    return_type: str = Field("entry", description="Return type for the base query (default 'entry').")

class RcsbPdbGetDateHistogramFacetsArgs(BaseModel):
    attribute_query_dict: Dict[str, Any] = Field(..., description="Dictionary for the base attribute query.")
    facet_name: str = Field(..., description="Descriptive name for the facet (e.g., 'Release Dates by Year').")
    facet_attribute: str = Field(..., description="Date attribute for histogram (e.g., 'rcsb_accession_info.initial_release_date').")
    interval: str = Field("year", description="Interval for date histogram (default 'year').")
    min_interval_population: int = Field(1, description="Minimum count for a bucket (default 1).")

class RcsbPdbGetCardinalityFacetArgs(BaseModel):
    attribute_query_dict: Dict[str, Any] = Field(..., description="Dictionary for the base attribute query.")
    facet_name: str = Field(..., description="Descriptive name for the facet (e.g., 'Unique Organism Count').")
    facet_attribute: str = Field(..., description="Attribute for which to count distinct values (e.g., 'rcsb_entity_source_organism.ncbi_scientific_name').")
    precision_threshold: int = Field(40000, description="Precision threshold for cardinality calculation (default 40000).")

class RcsbPdbCountQueryResultsArgs(BaseModel):
    query_type: str = Field(..., description="Type of query ('text', 'attribute', 'sequence', etc.).")
    query_params: Dict[str, Any] = Field(..., description="Parameters for the specified query type.")
    return_content_type: Optional[List[str]] = Field(None, description="Optional list to specify content types like [\"computational\", \"experimental\"].")

# --- RCSB PDB MCP Server Configuration ---
RCSB_MCP_HOST = os.getenv("RCSB_HOST", "localhost")
RCSB_MCP_PORT = os.getenv("RCSB_PORT", "8052")
RCSB_MCP_SERVER_URL = f"http://{RCSB_MCP_HOST}:{RCSB_MCP_PORT}/sse"
# rcsb_mcp_adapter = MCPAdapter(server_url=RCSB_MCP_SERVER_URL) # Removed
rcsb_router = APIRouter(prefix="/mcp/rcsb", tags=["RCSB PDB MCP Tools"])

# Explicit RCSB PDB tool endpoints
@rcsb_router.post("/text_search_pdb", summary="Text Search in RCSB PDB", description="Perform a text search in RCSB PDB.")
async def rcsb_text_search_pdb(tool_args: RcsbPdbTextSearchArgs = Body(...)):
    return await _run_mcp_tool_handler(RCSB_MCP_SERVER_URL, "text_search_pdb", tool_args.model_dump(exclude_none=True)) # Pass URL

@rcsb_router.post("/attribute_search_pdb", summary="Attribute Search in RCSB PDB", description="Perform an attribute search in RCSB PDB.")
async def rcsb_attribute_search_pdb(tool_args: RcsbPdbAttributeSearchArgs = Body(...)):
    return await _run_mcp_tool_handler(RCSB_MCP_SERVER_URL, "attribute_search_pdb", tool_args.model_dump(exclude_none=True)) # Pass URL

@rcsb_router.post("/combined_text_and_attribute_search", summary="Combined Text and Attribute Search", description="Combines a text query with multiple attribute queries.")
async def rcsb_combined_text_and_attribute_search(tool_args: RcsbPdbCombinedSearchArgs = Body(...)):
    return await _run_mcp_tool_handler(RCSB_MCP_SERVER_URL, "combined_text_and_attribute_search", tool_args.model_dump(exclude_none=True)) # Pass URL

@rcsb_router.post("/sequence_identity_search", summary="Sequence Identity Search", description="Find PDB entities with sequence similarity.")
async def rcsb_sequence_identity_search(tool_args: RcsbPdbSequenceIdentitySearchArgs = Body(...)):
    return await _run_mcp_tool_handler(RCSB_MCP_SERVER_URL, "sequence_identity_search", tool_args.model_dump(exclude_none=True)) # Pass URL

@rcsb_router.post("/sequence_motif_search", summary="Sequence Motif Search", description="Search for sequences containing a specific motif.")
async def rcsb_sequence_motif_search(tool_args: RcsbPdbSequenceMotifSearchArgs = Body(...)):
    return await _run_mcp_tool_handler(RCSB_MCP_SERVER_URL, "sequence_motif_search", tool_args.model_dump(exclude_none=True)) # Pass URL

@rcsb_router.post("/structure_similarity_by_entry_id", summary="Structure Similarity by Entry ID", description="Find structures similar to a given PDB entry.")
async def rcsb_structure_similarity_by_entry_id(tool_args: RcsbPdbStructSimilarityEntryIdArgs = Body(...)):
    return await _run_mcp_tool_handler(RCSB_MCP_SERVER_URL, "structure_similarity_by_entry_id", tool_args.model_dump(exclude_none=True)) # Pass URL

# @rcsb_router.post("/structure_similarity_by_file_url", summary="Structure Similarity by File URL", description="Find structures similar to a structure provided via a URL.")
# async def rcsb_structure_similarity_by_file_url(tool_args: RcsbPdbStructSimilarityFileUrlArgs = Body(...)):
#     return await _run_mcp_tool_handler(RCSB_MCP_SERVER_URL, "structure_similarity_by_file_url", tool_args.model_dump(exclude_none=True)) # Pass URL

@rcsb_router.post("/structure_motif_search_by_entry_id", summary="Structure Motif Search by Entry ID", description="Search for 3D structural motifs using a PDB entry as reference.")
async def rcsb_structure_motif_search_by_entry_id(tool_args: RcsbPdbStructMotifEntryIdArgs = Body(...)):
    return await _run_mcp_tool_handler(RCSB_MCP_SERVER_URL, "structure_motif_search_by_entry_id", tool_args.model_dump(exclude_none=True)) # Pass URL

# @rcsb_router.post("/get_term_facets", summary="Get Term Facets", description="Perform a query and get term-based aggregations (facets).")
# async def rcsb_get_term_facets(tool_args: RcsbPdbGetTermFacetsArgs = Body(...)):
#     return await _run_mcp_tool_handler(RCSB_MCP_SERVER_URL, "get_term_facets", tool_args.model_dump(exclude_none=True)) # Pass URL

# @rcsb_router.post("/get_histogram_facets", summary="Get Histogram Facets", description="Perform a query and get histogram-based aggregations for numeric attributes.")
# async def rcsb_get_histogram_facets(tool_args: RcsbPdbGetHistogramFacetsArgs = Body(...)):
#     return await _run_mcp_tool_handler(RCSB_MCP_SERVER_URL, "get_histogram_facets", tool_args.model_dump(exclude_none=True)) # Pass URL

# @rcsb_router.post("/get_date_histogram_facets", summary="Get Date Histogram Facets", description="Perform a query and get date histogram aggregations.")
# async def rcsb_get_date_histogram_facets(tool_args: RcsbPdbGetDateHistogramFacetsArgs = Body(...)):
#     return await _run_mcp_tool_handler(RCSB_MCP_SERVER_URL, "get_date_histogram_facets", tool_args.model_dump(exclude_none=True)) # Pass URL

# @rcsb_router.post("/get_cardinality_facet", summary="Get Cardinality Facet", description="Get the count of distinct values for a field.")
# async def rcsb_get_cardinality_facet(tool_args: RcsbPdbGetCardinalityFacetArgs = Body(...)):
#     return await _run_mcp_tool_handler(RCSB_MCP_SERVER_URL, "get_cardinality_facet", tool_args.model_dump(exclude_none=True)) # Pass URL

# @rcsb_router.post("/count_query_results", summary="Count Query Results", description="Counts the number of results for a given query.")
# async def rcsb_count_query_results(tool_args: RcsbPdbCountQueryResultsArgs = Body(...)):
#     return await _run_mcp_tool_handler(RCSB_MCP_SERVER_URL, "count_query_results", tool_args.model_dump(exclude_none=True)) # Pass URL

all_mcp_routers = [
    admet_router,
    alphafold_router,
    chembl_router,
    pubchem_router,
    rcsb_router,
    rxn_router,
    brics_router
]