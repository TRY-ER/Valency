/**
 * Redirection Mapper
 * Maps tool names to their corresponding redirection paths
 * This maps function call data.name values to specific application routes
 */

const redirectionMapper = {

    // admet tools
    "get_admet_prediction": "/optimization/admet",
    
    // alphafold tools
    "get_alphafold_prediction": "/explorer/alphafold",   
    "get_uniprot_summary": "/explorer/alphafold/summary",
    "get_alphafold_annotations": "/explorer/alphafold/annotations",

    // brics tools
    "get_brics_candidates": "/optimization/adv",

    // chembl tools 
    "find_molecule_by_pref_name": "/identification/chembl",
    "get_molecule_by_chembl_id": "/identification/chembl",
    "find_molecule_by_synonym": "/identification/chembl",
    "get_molecules_by_chembl_ids": "/identification/chembl",

    "find_similar_molecules_by_smiles": "/identification/chembl/similarity-search",
    "find_similar_molecules_by_chembl_id": "/identification/chembl/similarity-search",
    
    "get_activities_for_target": "/identification/chembl/activity-fetcher",
    "get_activities_for_molecule": "/identification/chembl/activity-fetcher",
    "find_target_by_gene_name": "/identification/chembl/utils", 
    "smiles_to_ctab": "/identification/chembl/utils",
    "compute_molecular_descriptors": "/identification/chembl/utils",
    "compute_structural_alerts": "/identification/chembl/utils",
    "standardize_molecule_from_smiles": "/identification/chembl/utils",
    "get_parent_molecule_from_smiles": "/identification/chembl/utils",


    // pubchem tools 
    "get_compound_by_cid": "/identification/pubchem",
    "get_cids_by_name": "/identification/pubchem",
    "get_cids_by_smiles": "/identification/pubchem",
    "get_cids_by_inchikey": "/identification/pubchem",
    "get_compound_properties": "/identification/pubchem/utils",
    "get_compound_synonyms_by_cid": "/identification/pubchem/utils",
    "fast_identity_search_by_cid": "/identification/pubchem/utils",
    "fast_similarity_2d_search_by_cid": "/identification/pubchem/similarity-search",
    "fast_substructure_search_by_smiles": "/identification/pubchem/utils", 


    // rcsb tools
    "text_search_pdb": "/explorer/rcsb",
    "get_protein_details_by_id_pypdb": "/explorer/rcsb",
    "structure_similarity_by_entry_id":  "/explorer/rcsb/ssearch",

    // uniprot tools
    "search_uniprotkb": "/explorer/uniprot",
    "get_uniprotkb_entry": "/explorer/uniprot", 

    // Add more tool mappings here as needed
    // "another_tool": "/another/route",
};

/**
 * Gets the redirection path for a given tool name
 * @param {string} toolName - The name of the tool
 * @returns {string|null} - The redirection path or null if no mapping exists
 */
export const getRedirectionPath = (toolName) => {
    return redirectionMapper[toolName] || null;
};

/**
 * Checks if a tool has a redirection path configured
 * @param {string} toolName - The name of the tool
 * @returns {boolean} - True if redirection is available, false otherwise
 */
export const hasRedirection = (toolName) => {
    return toolName in redirectionMapper;
};

/**
 * Gets all available redirections (for debugging/admin purposes)
 * @returns {Object} - The complete redirection mapping
 */
export const getAllRedirections = () => {
    return redirectionMapper;
};

export default {
    getRedirectionPath,
    hasRedirection,
    getAllRedirections
};
