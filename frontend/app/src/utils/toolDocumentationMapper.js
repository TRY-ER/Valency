/**
 * Tool Documentation Mapper
 * Maps tool names to their corresponding markdown documentation file paths
 */

const toolDocumentationMapper = {
    // Molecular tools
    "get_admet_prediction": "/markdown_repo/ADMETPredictor.md",
    
    // AlphaFold tools
    "get_alphafold_prediction": "/markdown_repo/UniProtViewer.md",
    "get_uniprot_summary": "/markdown_repo/UniProtSummaryViewer.md",
    "get_alphafold_annotations": "/markdown_repo/AlphaFoldAnnotationsViewer.md",
    
    // BRICS tools
    "get_brics_candidates": "/markdown_repo/AdvancedBRICS.md",
    
    // ChEMBL tools
    "find_molecule_by_pref_name": "/markdown_repo/ChEMBLGetter.md",
    "get_molecule_by_chembl_id": "/markdown_repo/ChEMBLGetter.md",
    "find_molecule_by_synonym": "/markdown_repo/ChEMBLGetter.md",
    "get_molecules_by_chembl_ids": "/markdown_repo/ChEMBLGetter.md",
    "find_similar_molecules_by_smiles": "/markdown_repo/ChEMBLSimilarityGetter.md",
    "find_similar_molecules_by_chembl_id": "/markdown_repo/ChEMBLSimilarityGetter.md",
    "get_activities_for_target": "/markdown_repo/ChEMBLActivityFetcher.md",
    "get_activities_for_molecule": "/markdown_repo/ChEMBLActivityFetcher.md",
    
    // ChEMBL Utilities
    "find_target_by_gene_name": "/markdown_repo/ChEMBLUtilities.md",
    "smiles_to_ctab": "/markdown_repo/ChEMBLUtilities.md",
    "compute_molecular_descriptors": "/markdown_repo/ChEMBLUtilities.md",
    "compute_structural_alerts": "/markdown_repo/ChEMBLUtilities.md",
    "standardize_molecule_from_smiles": "/markdown_repo/ChEMBLUtilities.md",
    "get_parent_molecule_from_smiles": "/markdown_repo/ChEMBLUtilities.md",
    
    // PubChem tools
    "get_compound_by_cid": "/markdown_repo/PubChemGetter.md",
    "get_cids_by_name": "/markdown_repo/PubChemGetter.md",
    "get_cids_by_smiles": "/markdown_repo/PubChemGetter.md",
    "get_cids_by_inchikey": "/markdown_repo/PubChemGetter.md",
    "get_compound_properties": "/markdown_repo/PubChemUtilities.md",
    "get_compound_synonyms_by_cid": "/markdown_repo/PubChemUtilities.md",
    "fast_identity_search_by_cid": "/markdown_repo/PubChemUtilities.md",
    "fast_similarity_2d_search_by_cid": "/markdown_repo/PubChemSimilarityGetter.md",
    "fast_substructure_search_by_smiles": "/markdown_repo/PubChemUtilities.md",
    
    // RCSB PDB tools
    "text_search_pdb": "/markdown_repo/RCSBPDBExplorer.md",
    "get_protein_details_by_id_pypdb": "/markdown_repo/RCSBPDBExplorer.md",
    "structure_similarity_by_entry_id": "/markdown_repo/RCSBStructureSimilaritySearch.md",
    
    // UniProt tools
    "search_uniprotkb": "/markdown_repo/UniProtExplorer.md",
    "get_uniprotkb_entry": "/markdown_repo/UniProtExplorer.md",
    
    // LSTM tools
    "lstm_psmiles_generation": "/markdown_repo/LSTMPSMILES.md",
    "lstm_wdg_generation": "/markdown_repo/LSTMWDG.md",
    
    // BRICS generators
    "brics_smiles_generation": "/markdown_repo/BRICSSMILES.md",
    "brics_psmiles_generation": "/markdown_repo/BRICSPSMILES.md",
    
    // Explorer tools
    "molecule_explorer": "/markdown_repo/MoleculeExplorer.md",
    "protein_explorer": "/markdown_repo/ProteinExplorer.md",
    "polymer_explorer": "/markdown_repo/PolymerExplorer.md",
    
    // Similarity search tools
    "molecule_similarity_search": "/markdown_repo/MoleculeSimilaritySearch.md",
    "polymer_similarity_search": "/markdown_repo/PolymerSimilaritySearch.md",
    "protein_similarity_search": "/markdown_repo/ProteinSimilaritySearch.md",
};

/**
 * Get documentation file path for a tool
 * @param {string} toolName - The name of the tool
 * @returns {string|null} - The documentation file path or null if not found
 */
export const getToolDocumentationPath = (toolName) => {
    return toolDocumentationMapper[toolName] || null;
};

/**
 * Check if a tool has documentation available
 * @param {string} toolName - The name of the tool
 * @returns {boolean} - Whether documentation is available
 */
export const hasDocumentation = (toolName) => {
    return toolName in toolDocumentationMapper;
};

/**
 * Get all tools that have documentation
 * @returns {Array<string>} - Array of tool names with documentation
 */
export const getToolsWithDocumentation = () => {
    return Object.keys(toolDocumentationMapper);
};

export default toolDocumentationMapper;
