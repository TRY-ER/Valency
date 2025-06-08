/**
 * Tool Type Mapper
 * Maps tool names to their corresponding type annotations
 * If a tool name has null value, no type is added
 * If a tool name has a string value, that type is added to the data structure
 */

const toolTypeMapper = {
    "get_admet_prediction": null,
    "get_alphafold_prediction": null,   
    "get_uniprot_summary": null,
    "get_alphafold_annotations": null,
    "get_brics_candidates": null,
    "find_molecule_by_pref_name" : "pref_name",
    "get_molecule_by_chembl_id" : "chembl_id",
    "find_molecule_by_synonym" : "synonym",
    "get_molecules_by_chembl_ids": "multiple_chembl_ids",
    "find_similar_molecules_by_smiles" : "smiles",
    "find_similar_molecules_by_chembl_id" : "chembl_id",
    "get_activities_for_target": "target",
    "get_activities_for_molecule": "molecule",
    "find_target_by_gene_name": null,
    "smiles_to_ctab": null,
    "compute_molecular_descriptors": null,
    "compute_structural_alerts": null,
    "standardize_molecule_from_smiles": null,
    "get_parent_molecule_from_smiles": null,
    "get_compound_by_cid": "cid",
    "get_cids_by_name": "name",
    "get_cids_by_smiles": "smiles",
    "get_cids_by_inchikey": "inchikey",
    "get_compound_properties": null,
    "get_compound_synonyms_by_cid": null,
    "fast_identity_search_by_cid": "cid",
    "fast_similarity_2d_search_by_cid": "cid",
    "fast_substructure_search_by_smiles": null
    // Add more tool mappings here as needed
    // "another_tool": "another_type",
    // "simple_tool": null,
};

/**
 * Processes tool data based on tool name and adds type annotation if configured
 * @param {string} toolName - The name of the tool
 * @param {any} toolData - The raw tool data
 * @returns {any} - Processed tool data with or without type annotation
 */
const processToolData = (toolName, toolData) => {
    // Check if tool name exists in mapper
    if (toolName in toolTypeMapper) {
        const typeAnnotation = toolTypeMapper[toolName];
        
        // If type annotation is null, return data as-is
        if (typeAnnotation === null) {
            return toolData;
        }
        
        // If type annotation exists, wrap data with type
        return {
            type: typeAnnotation,
            content: toolData
        };
    }
    
    // If tool name not in mapper, return data as-is
    return toolData;
};

/**
 * Adds a new tool type mapping
 * @param {string} toolName - The name of the tool
 * @param {string|null} typeAnnotation - The type annotation or null
 */
export const addToolTypeMapping = (toolName, typeAnnotation) => {
    toolTypeMapper[toolName] = typeAnnotation;
};

/**
 * Gets the current tool type mapper for debugging/inspection
 * @returns {object} - The current tool type mapper
 */
export const getToolTypeMapper = () => {
    return { ...toolTypeMapper };
};

export default { processToolData, addToolTypeMapping, getToolTypeMapper };
