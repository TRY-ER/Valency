import axios from 'axios';

const MCP_TOOLS_BASE_URL = 'http://localhost:8000'; // Base URL for MCP tools

/**
 * Retrieves the auth token from localStorage.
 * In a real app, this might come from a more robust auth service or context.
 */
const getAuthToken = () => {
  return localStorage.getItem('access_token');
};

/**
 * Creates an axios instance with common configurations for the MCP tools service.
 */
const mcpToolsApiClient = axios.create({
  baseURL: MCP_TOOLS_BASE_URL,
});

/**
 * Interceptor to add the Authorization header to requests if a token exists.
 */
mcpToolsApiClient.interceptors.request.use(
  (config) => {
    const token = getAuthToken();
    if (token) {
      config.headers['Authorization'] = `Bearer ${token}`;
    }
    return config;
  },
  (error) => {
    return Promise.reject(error);
  }
);

/**
 * Generic function to call MCP tools endpoints.
 * @param {string} toolPath - The tool path (e.g., '/mcp/admet/get_admet_prediction')
 * @param {object} toolArgs - The arguments to pass to the tool
 * @param {string} method - HTTP method (default: 'POST')
 * @returns {Promise<object>} The response data from the tool endpoint
 */
const callMcpTool = async (toolPath, toolArgs = {}, method = 'POST') => {
  try {
    let response;
    
    if (method.toUpperCase() === 'GET') {
      response = await mcpToolsApiClient.get(toolPath, { params: toolArgs });
    } else if (method.toUpperCase() === 'POST') {
      response = await mcpToolsApiClient.post(toolPath, toolArgs);
    } else if (method.toUpperCase() === 'PUT') {
      response = await mcpToolsApiClient.put(toolPath, toolArgs);
    } else if (method.toUpperCase() === 'DELETE') {
      response = await mcpToolsApiClient.delete(toolPath, { data: toolArgs });
    } else {
      throw new Error(`Unsupported HTTP method: ${method}`);
    }
    
    return response.data;
  } catch (error) {
    console.error(`Error calling MCP tool ${toolPath}:`, error.response ? error.response.data : error.message);
    
    // Handle specific error cases
    if (error.response) {
      const { status, data } = error.response;
      
      if (status === 401) {
        console.error('Authentication failed for MCP tool call');
        // You might want to trigger a re-login here
      } else if (status === 404) {
        console.error(`MCP tool endpoint not found: ${toolPath}`);
      } else if (status === 422) {
        console.error('Invalid arguments provided to MCP tool:', data);
      } else if (status >= 500) {
        console.error('Server error occurred while calling MCP tool');
      }
    }
    
    throw error;
  }
};

// ==================== ADMET TOOLS ====================

/**
 * @typedef {Object} AdmetPredictionArgs
 * @property {string} smiles - The SMILES string representing the molecule.
 */

/**
 * Get ADMET (Absorption, Distribution, Metabolism, Excretion, Toxicity) predictions for a given SMILES string.
 * @param {AdmetPredictionArgs} args - The arguments for ADMET prediction
 * @returns {Promise<object>} A promise that resolves to the ADMET prediction data
 */
export const getAdmetPrediction = async (args) => {
  if (!args || !args.smiles) {
    console.error('getAdmetPrediction: smiles is required');
    return Promise.reject(new Error('SMILES string is required'));
  }
  
  if (typeof args.smiles !== 'string' || args.smiles.trim() === '') {
    console.error('getAdmetPrediction: smiles must be a non-empty string');
    return Promise.reject(new Error('SMILES must be a non-empty string'));
  }
  
  try {
    const result = await callMcpTool('/mcp/admet/get_admet_prediction', args);
    console.log('ADMET prediction retrieved successfully');
    return result;
  } catch (error) {
    console.error('Failed to get ADMET prediction:', error.message);
    throw error;
  }
};

// ==================== ALPHAFOLD TOOLS ====================

/**
 * @typedef {Object} AlphafoldPredictionArgs
 * @property {string} uniprot_accession - The UniProt accession ID for the protein.
 */

/**
 * Get AlphaFold prediction models for a UniProt accession.
 * @param {AlphafoldPredictionArgs} args - The arguments for AlphaFold prediction
 * @returns {Promise<object>} A promise that resolves to the AlphaFold prediction data
 */
export const getAlphafoldPrediction = async (args) => {
  if (!args || !args.uniprot_accession) {
    console.error('getAlphafoldPrediction: uniprot_accession is required');
    return Promise.reject(new Error('UniProt accession is required'));
  }
  
  if (typeof args.uniprot_accession !== 'string' || args.uniprot_accession.trim() === '') {
    console.error('getAlphafoldPrediction: uniprot_accession must be a non-empty string');
    return Promise.reject(new Error('UniProt accession must be a non-empty string'));
  }
  
  try {
    const result = await callMcpTool('/mcp/alphafold/get_alphafold_prediction', args);
    console.log('AlphaFold prediction retrieved successfully');
    return result;
  } catch (error) {
    console.error('Failed to get AlphaFold prediction:', error.message);
    throw error;
  }
};

/**
 * @typedef {Object} UniprotSummaryArgs
 * @property {string} qualifier - UniProtKB accession number (AC), entry name (ID), or CRC64 checksum of the UniProt sequence (e.g., 'Q5VSL9').
 * @property {number} [start] - The start position of the residue range (optional).
 * @property {number} [end] - The end position of the residue range (optional).
 */

/**
 * Get UniProt summary details for a residue range.
 * @param {UniprotSummaryArgs} args - The arguments for UniProt summary
 * @returns {Promise<object>} A promise that resolves to the UniProt summary data
 */
export const getUniprotSummary = async (args) => {
  if (!args || !args.qualifier) {
    console.error('getUniprotSummary: qualifier is required');
    return Promise.reject(new Error('UniProt qualifier is required'));
  }
  
  if (typeof args.qualifier !== 'string' || args.qualifier.trim() === '') {
    console.error('getUniprotSummary: qualifier must be a non-empty string');
    return Promise.reject(new Error('UniProt qualifier must be a non-empty string'));
  }
  
  try {
    const result = await callMcpTool('/mcp/alphafold/get_uniprot_summary', args);
    console.log('UniProt summary retrieved successfully');
    return result;
  } catch (error) {
    console.error('Failed to get UniProt summary:', error.message);
    throw error;
  }
};

/**
 * @typedef {Object} AlphafoldAnnotationsArgs
 * @property {string} qualifier - UniProt accession (e.g., 'Q5VSL9').
 * @property {string} annotation_type - Type of annotation (e.g., 'MUTAGEN' for AlphaMissense).
 */

/**
 * Get AlphaFold annotations for a specific UniProt accession and annotation type.
 * @param {AlphafoldAnnotationsArgs} args - The arguments for AlphaFold annotations
 * @returns {Promise<object>} A promise that resolves to the AlphaFold annotations data
 */
export const getAlphafoldAnnotations = async (args) => {
  if (!args || !args.qualifier) {
    console.error('getAlphafoldAnnotations: qualifier is required');
    return Promise.reject(new Error('UniProt qualifier is required'));
  }
  
  if (!args.annotation_type) {
    console.error('getAlphafoldAnnotations: annotation_type is required');
    return Promise.reject(new Error('Annotation type is required'));
  }
  
  if (typeof args.qualifier !== 'string' || args.qualifier.trim() === '') {
    console.error('getAlphafoldAnnotations: qualifier must be a non-empty string');
    return Promise.reject(new Error('UniProt qualifier must be a non-empty string'));
  }
  
  if (typeof args.annotation_type !== 'string' || args.annotation_type.trim() === '') {
    console.error('getAlphafoldAnnotations: annotation_type must be a non-empty string');
    return Promise.reject(new Error('Annotation type must be a non-empty string'));
  }
  
  try {
    const result = await callMcpTool('/mcp/alphafold/get_alphafold_annotations', args);
    console.log('AlphaFold annotations retrieved successfully');
    return result;
  } catch (error) {
    console.error('Failed to get AlphaFold annotations:', error.message);
    throw error;
  }
};

// ==================== BRICS TOOLS ====================

/**
 * @typedef {Object} BricsCandidatesArgs
 * @property {string[]} smiles_list - A list of SMILES strings.
 * @property {boolean} [is_polymer=false] - A boolean flag indicating if the input SMILES are polymers.
 */

/**
 * Generate molecular candidates from a list of SMILES strings using the BRICS algorithm.
 * @param {BricsCandidatesArgs} args - The arguments for BRICS candidates generation
 * @returns {Promise<object>} A promise that resolves to the BRICS candidates data
 */
export const getBricsCandidates = async (args) => {
  if (!args || !args.smiles_list) {
    console.error('getBricsCandidates: smiles_list is required');
    return Promise.reject(new Error('SMILES list is required'));
  }
  
  if (!Array.isArray(args.smiles_list) || args.smiles_list.length === 0) {
    console.error('getBricsCandidates: smiles_list must be a non-empty array');
    return Promise.reject(new Error('SMILES list must be a non-empty array'));
  }
  
  // Validate that all elements in the array are strings
  if (!args.smiles_list.every(smiles => typeof smiles === 'string' && smiles.trim() !== '')) {
    console.error('getBricsCandidates: all SMILES strings must be non-empty strings');
    return Promise.reject(new Error('All SMILES strings must be non-empty strings'));
  }
  
  try {
    const result = await callMcpTool('/mcp/brics/get_brics_candidates', args);
    console.log('BRICS candidates retrieved successfully');
    return result;
  } catch (error) {
    console.error('Failed to get BRICS candidates:', error.message);
    throw error;
  }
};

// ==================== CHEMBL TOOLS ====================

/**
 * @typedef {Object} ChemblGetMoleculeByIdArgs
 * @property {string} chembl_id - The ChEMBL ID of the molecule (e.g., 'CHEMBL192').
 */

/**
 * Retrieve a specific molecule by its unique ChEMBL ID.
 * @param {ChemblGetMoleculeByIdArgs} args - The arguments for molecule retrieval
 * @returns {Promise<object>} A promise that resolves to the molecule data
 */
export const getMoleculeByChemblId = async (args) => {
  if (!args || !args.chembl_id) {
    console.error('getMoleculeByChemblId: chembl_id is required');
    return Promise.reject(new Error('ChEMBL ID is required'));
  }
  
  if (typeof args.chembl_id !== 'string' || args.chembl_id.trim() === '') {
    console.error('getMoleculeByChemblId: chembl_id must be a non-empty string');
    return Promise.reject(new Error('ChEMBL ID must be a non-empty string'));
  }
  
  try {
    const result = await callMcpTool('/mcp/chembl/get_molecule_by_chembl_id', args);
    console.log('ChEMBL molecule retrieved successfully');
    return result;
  } catch (error) {
    console.error('Failed to get ChEMBL molecule:', error.message);
    throw error;
  }
};

/**
 * @typedef {Object} ChemblFindMoleculeByPrefNameArgs
 * @property {string} pref_name - The preferred name to search for (e.g., 'Aspirin').
 */

/**
 * Search for molecules by their preferred name.
 * @param {ChemblFindMoleculeByPrefNameArgs} args - The arguments for molecule search
 * @returns {Promise<object>} A promise that resolves to the molecule search results
 */
export const findMoleculeByPrefName = async (args) => {
  if (!args || !args.pref_name) {
    console.error('findMoleculeByPrefName: pref_name is required');
    return Promise.reject(new Error('Preferred name is required'));
  }
  
  if (typeof args.pref_name !== 'string' || args.pref_name.trim() === '') {
    console.error('findMoleculeByPrefName: pref_name must be a non-empty string');
    return Promise.reject(new Error('Preferred name must be a non-empty string'));
  }
  
  try {
    const result = await callMcpTool('/mcp/chembl/find_molecule_by_pref_name', args);
    console.log('ChEMBL molecule search by name completed successfully');
    return result;
  } catch (error) {
    console.error('Failed to find ChEMBL molecule by name:', error.message);
    throw error;
  }
};

/**
 * @typedef {Object} ChemblFindMoleculeBySynonymArgs
 * @property {string} synonym - The synonym to search for (case-insensitive exact match).
 */

/**
 * Find molecules by their synonyms.
 * @param {ChemblFindMoleculeBySynonymArgs} args - The arguments for synonym search
 * @returns {Promise<object>} A promise that resolves to the molecule search results
 */
export const findMoleculeBySynonym = async (args) => {
  if (!args || !args.synonym) {
    console.error('findMoleculeBySynonym: synonym is required');
    return Promise.reject(new Error('Synonym is required'));
  }
  
  if (typeof args.synonym !== 'string' || args.synonym.trim() === '') {
    console.error('findMoleculeBySynonym: synonym must be a non-empty string');
    return Promise.reject(new Error('Synonym must be a non-empty string'));
  }
  
  try {
    const result = await callMcpTool('/mcp/chembl/find_molecule_by_synonym', args);
    console.log('ChEMBL molecule search by synonym completed successfully');
    return result;
  } catch (error) {
    console.error('Failed to find ChEMBL molecule by synonym:', error.message);
    throw error;
  }
};

/**
 * @typedef {Object} ChemblGetMoleculesByIdsArgs
 * @property {string[]} chembl_ids - A list of ChEMBL IDs (e.g., ['CHEMBL25', 'CHEMBL192']).
 */

/**
 * Retrieve multiple molecules by a list of ChEMBL IDs.
 * @param {ChemblGetMoleculesByIdsArgs} args - The arguments for multiple molecule retrieval
 * @returns {Promise<object>} A promise that resolves to the molecules data
 */
export const getMoleculesByChemblIds = async (args) => {
  if (!args || !args.chembl_ids) {
    console.error('getMoleculesByChemblIds: chembl_ids is required');
    return Promise.reject(new Error('ChEMBL IDs list is required'));
  }
  
  if (!Array.isArray(args.chembl_ids) || args.chembl_ids.length === 0) {
    console.error('getMoleculesByChemblIds: chembl_ids must be a non-empty array');
    return Promise.reject(new Error('ChEMBL IDs must be a non-empty array'));
  }
  
  if (!args.chembl_ids.every(id => typeof id === 'string' && id.trim() !== '')) {
    console.error('getMoleculesByChemblIds: all ChEMBL IDs must be non-empty strings');
    return Promise.reject(new Error('All ChEMBL IDs must be non-empty strings'));
  }
  
  try {
    const result = await callMcpTool('/mcp/chembl/get_molecules_by_chembl_ids', args);
    console.log('ChEMBL molecules retrieved successfully');
    return result;
  } catch (error) {
    console.error('Failed to get ChEMBL molecules:', error.message);
    throw error;
  }
};

/**
 * @typedef {Object} ChemblGetMoleculeImageSvgArgs
 * @property {string} chembl_id - The ChEMBL ID of the molecule.
 */

/**
 * Get the 2D chemical structure image of a molecule as an SVG string.
 * @param {ChemblGetMoleculeImageSvgArgs} args - The arguments for molecule image retrieval
 * @returns {Promise<object>} A promise that resolves to the SVG image data
 */
export const getMoleculeImageSvg = async (args) => {
  if (!args || !args.chembl_id) {
    console.error('getMoleculeImageSvg: chembl_id is required');
    return Promise.reject(new Error('ChEMBL ID is required'));
  }
  
  if (typeof args.chembl_id !== 'string' || args.chembl_id.trim() === '') {
    console.error('getMoleculeImageSvg: chembl_id must be a non-empty string');
    return Promise.reject(new Error('ChEMBL ID must be a non-empty string'));
  }
  
  try {
    const result = await callMcpTool('/mcp/chembl/get_molecule_image_svg', args);
    console.log('ChEMBL molecule image retrieved successfully');
    return result;
  } catch (error) {
    console.error('Failed to get ChEMBL molecule image:', error.message);
    throw error;
  }
};

/**
 * @typedef {Object} ChemblFindSimilarMoleculesBySmilesArgs
 * @property {string} smiles - The SMILES string of the query molecule (e.g., 'CCO').
 * @property {number} [similarity_threshold=70] - Minimum similarity percentage (0-100).
 */

/**
 * Find molecules structurally similar to a given SMILES string.
 * @param {ChemblFindSimilarMoleculesBySmilesArgs} args - The arguments for similarity search
 * @returns {Promise<object>} A promise that resolves to the similar molecules data
 */
export const findSimilarMoleculesBySmiles = async (args) => {
  if (!args || !args.smiles) {
    console.error('findSimilarMoleculesBySmiles: smiles is required');
    return Promise.reject(new Error('SMILES string is required'));
  }
  
  if (typeof args.smiles !== 'string' || args.smiles.trim() === '') {
    console.error('findSimilarMoleculesBySmiles: smiles must be a non-empty string');
    return Promise.reject(new Error('SMILES string must be a non-empty string'));
  }
  
  try {
    const result = await callMcpTool('/mcp/chembl/find_similar_molecules_by_smiles', args);
    console.log('ChEMBL similar molecules search completed successfully');
    return result;
  } catch (error) {
    console.error('Failed to find similar ChEMBL molecules by SMILES:', error.message);
    throw error;
  }
};

/**
 * @typedef {Object} ChemblFindSimilarMoleculesByIdArgs
 * @property {string} chembl_id - The ChEMBL ID of the query molecule.
 * @property {number} [similarity_threshold=70] - Minimum similarity percentage (0-100).
 */

/**
 * Find molecules structurally similar to a given ChEMBL ID.
 * @param {ChemblFindSimilarMoleculesByIdArgs} args - The arguments for similarity search
 * @returns {Promise<object>} A promise that resolves to the similar molecules data
 */
export const findSimilarMoleculesByChemblId = async (args) => {
  if (!args || !args.chembl_id) {
    console.error('findSimilarMoleculesByChemblId: chembl_id is required');
    return Promise.reject(new Error('ChEMBL ID is required'));
  }
  
  if (typeof args.chembl_id !== 'string' || args.chembl_id.trim() === '') {
    console.error('findSimilarMoleculesByChemblId: chembl_id must be a non-empty string');
    return Promise.reject(new Error('ChEMBL ID must be a non-empty string'));
  }
  
  try {
    const result = await callMcpTool('/mcp/chembl/find_similar_molecules_by_chembl_id', args);
    console.log('ChEMBL similar molecules search completed successfully');
    return result;
  } catch (error) {
    console.error('Failed to find similar ChEMBL molecules by ID:', error.message);
    throw error;
  }
};

/**
 * @typedef {Object} ChemblGetApprovedDrugsArgs
 * @property {boolean} [order_by_mw=false] - If true, sorts results by molecular weight.
 */

/**
 * Retrieve all drugs that have reached the maximum clinical trial phase.
 * @param {ChemblGetApprovedDrugsArgs} args - The arguments for approved drugs retrieval
 * @returns {Promise<object>} A promise that resolves to the approved drugs data
 */
export const getApprovedDrugs = async (args = {}) => {
  try {
    const result = await callMcpTool('/mcp/chembl/get_approved_drugs', args);
    console.log('ChEMBL approved drugs retrieved successfully');
    return result;
  } catch (error) {
    console.error('Failed to get ChEMBL approved drugs:', error.message);
    throw error;
  }
};

/**
 * @typedef {Object} ChemblGetActivitiesForTargetArgs
 * @property {string} biological_target_name - Biological target name (e.g., 'hERG').
 * @property {string} [standard_type='IC50'] - Bioactivity measurement type (e.g., 'IC50', 'Ki').
 */

/**
 * Fetch bioactivity data for a specific biological target.
 * @param {ChemblGetActivitiesForTargetArgs} args - The arguments for target activities retrieval
 * @returns {Promise<object>} A promise that resolves to the target activities data
 */
export const getActivitiesForTarget = async (args) => {
  if (!args || !args.biological_target_name) {
    console.error('getActivitiesForTarget: target_chembl_id is required');
    return Promise.reject(new Error('Target ChEMBL ID is required'));
  }
  
  if (typeof args.biological_target_name!== 'string' || args.biological_target_name.trim() === '') {
    console.error('getActivitiesForTarget: target_chembl_id must be a non-empty string');
    return Promise.reject(new Error('Target ChEMBL ID must be a non-empty string'));
  }
  
  try {
    const result = await callMcpTool('/mcp/chembl/get_activities_for_target', args);
    console.log('ChEMBL target activities retrieved successfully');
    return result;
  } catch (error) {
    console.error('Failed to get ChEMBL target activities:', error.message);
    throw error;
  }
};

/**
 * @typedef {Object} ChemblGetActivitiesForMoleculeArgs
 * @property {string} molecule_chembl_id - ChEMBL ID of the molecule.
 * @property {boolean} [pchembl_value_exists=true] - If true, only returns activities with a pChEMBL value.
 */

/**
 * Retrieve all recorded bioactivities for a specific molecule.
 * @param {ChemblGetActivitiesForMoleculeArgs} args - The arguments for molecule activities retrieval
 * @returns {Promise<object>} A promise that resolves to the molecule activities data
 */
export const getActivitiesForMolecule = async (args) => {
  if (!args || !args.molecule_chembl_id) {
    console.error('getActivitiesForMolecule: molecule_chembl_id is required');
    return Promise.reject(new Error('Molecule ChEMBL ID is required'));
  }
  
  if (typeof args.molecule_chembl_id !== 'string' || args.molecule_chembl_id.trim() === '') {
    console.error('getActivitiesForMolecule: molecule_chembl_id must be a non-empty string');
    return Promise.reject(new Error('Molecule ChEMBL ID must be a non-empty string'));
  }
  
  try {
    const result = await callMcpTool('/mcp/chembl/get_activities_for_molecule', args);
    console.log('ChEMBL molecule activities retrieved successfully');
    return result;
  } catch (error) {
    console.error('Failed to get ChEMBL molecule activities:', error.message);
    throw error;
  }
};

/**
 * @typedef {Object} ChemblFindTargetByGeneNameArgs
 * @property {string} gene_name - Gene name or symbol (e.g., 'BRCA1'). Case-insensitive contains match.
 */

/**
 * Search for biological targets by a gene name or symbol.
 * @param {ChemblFindTargetByGeneNameArgs} args - The arguments for target search by gene name
 * @returns {Promise<object>} A promise that resolves to the target search results
 */
export const findTargetByGeneName = async (args) => {
  if (!args || !args.gene_name) {
    console.error('findTargetByGeneName: gene_name is required');
    return Promise.reject(new Error('Gene name is required'));
  }
  
  if (typeof args.gene_name !== 'string' || args.gene_name.trim() === '') {
    console.error('findTargetByGeneName: gene_name must be a non-empty string');
    return Promise.reject(new Error('Gene name must be a non-empty string'));
  }
  
  try {
    const result = await callMcpTool('/mcp/chembl/find_target_by_gene_name', args);
    console.log('ChEMBL target search by gene name completed successfully');
    return result;
  } catch (error) {
    console.error('Failed to find ChEMBL target by gene name:', error.message);
    throw error;
  }
};

/**
 * @typedef {Object} ChemblGetMoleculesByFilterArgs
 * @property {Object} filters - Django-style filters (e.g., {'molecule_properties__mw_freebase__gte': 200}).
 * @property {string[]} [only_fields] - Optional list of specific fields to return.
 * @property {string[]} [order_by] - Fields to sort by (e.g., ['molecule_properties__mw_freebase']).
 */

/**
 * Retrieve molecules based on a custom set of filter conditions.
 * @param {ChemblGetMoleculesByFilterArgs} args - The arguments for filtered molecule retrieval
 * @returns {Promise<object>} A promise that resolves to the filtered molecules data
 */
export const getMoleculesByFilter = async (args) => {
  if (!args || !args.filters) {
    console.error('getMoleculesByFilter: filters is required');
    return Promise.reject(new Error('Filters object is required'));
  }
  
  if (typeof args.filters !== 'object' || args.filters === null) {
    console.error('getMoleculesByFilter: filters must be an object');
    return Promise.reject(new Error('Filters must be an object'));
  }
  
  try {
    const result = await callMcpTool('/mcp/chembl/get_molecules_by_filter', args);
    console.log('ChEMBL filtered molecules retrieved successfully');
    return result;
  } catch (error) {
    console.error('Failed to get ChEMBL molecules by filter:', error.message);
    throw error;
  }
};

/**
 * @typedef {Object} ChemblGetActivitiesByFilterArgs
 * @property {Object} filters - Filters for activity fields (e.g., {'pchembl_value__gte': 6.0}).
 * @property {string[]} [only_fields] - Optional list of specific fields to return.
 * @property {string[]} [order_by] - Fields to sort by (e.g., ['-pchembl_value']).
 */

/**
 * Retrieve bioactivity data based on a custom set of filter conditions.
 * @param {ChemblGetActivitiesByFilterArgs} args - The arguments for filtered activities retrieval
 * @returns {Promise<object>} A promise that resolves to the filtered activities data
 */
export const getActivitiesByFilter = async (args) => {
  if (!args || !args.filters) {
    console.error('getActivitiesByFilter: filters is required');
    return Promise.reject(new Error('Filters object is required'));
  }
  
  if (typeof args.filters !== 'object' || args.filters === null) {
    console.error('getActivitiesByFilter: filters must be an object');
    return Promise.reject(new Error('Filters must be an object'));
  }
  
  try {
    const result = await callMcpTool('/mcp/chembl/get_activities_by_filter', args);
    console.log('ChEMBL filtered activities retrieved successfully');
    return result;
  } catch (error) {
    console.error('Failed to get ChEMBL activities by filter:', error.message);
    throw error;
  }
};

/**
 * @typedef {Object} ChemblGetTargetsByFilterArgs
 * @property {Object} filters - Filters for target fields (e.g., {'target_type': 'SINGLE PROTEIN'}).
 * @property {string[]} [only_fields] - Optional list of specific fields to return.
 * @property {string[]} [order_by] - Fields to sort by (e.g., ['pref_name']).
 */

/**
 * Retrieve biological targets based on a custom set of filter conditions.
 * @param {ChemblGetTargetsByFilterArgs} args - The arguments for filtered targets retrieval
 * @returns {Promise<object>} A promise that resolves to the filtered targets data
 */
export const getTargetsByFilter = async (args) => {
  if (!args || !args.filters) {
    console.error('getTargetsByFilter: filters is required');
    return Promise.reject(new Error('Filters object is required'));
  }
  
  if (typeof args.filters !== 'object' || args.filters === null) {
    console.error('getTargetsByFilter: filters must be an object');
    return Promise.reject(new Error('Filters must be an object'));
  }
  
  try {
    const result = await callMcpTool('/mcp/chembl/get_targets_by_filter', args);
    console.log('ChEMBL filtered targets retrieved successfully');
    return result;
  } catch (error) {
    console.error('Failed to get ChEMBL targets by filter:', error.message);
    throw error;
  }
};

/**
 * @typedef {Object} ChemblSmilesToCtabArgs
 * @property {string} smiles - SMILES string to convert.
 */

/**
 * Convert a SMILES string to a CTAB block.
 * @param {ChemblSmilesToCtabArgs} args - The arguments for SMILES to CTAB conversion
 * @returns {Promise<object>} A promise that resolves to the CTAB data
 */
export const smilesToCtab = async (args) => {
  if (!args || !args.smiles) {
    console.error('smilesToCtab: smiles is required');
    return Promise.reject(new Error('SMILES string is required'));
  }
  
  if (typeof args.smiles !== 'string' || args.smiles.trim() === '') {
    console.error('smilesToCtab: smiles must be a non-empty string');
    return Promise.reject(new Error('SMILES string must be a non-empty string'));
  }
  
  try {
    const result = await callMcpTool('/mcp/chembl/smiles_to_ctab', args);
    console.log('ChEMBL SMILES to CTAB conversion completed successfully');
    return result;
  } catch (error) {
    console.error('Failed to convert SMILES to CTAB:', error.message);
    throw error;
  }
};

/**
 * @typedef {Object} ChemblComputeMolecularDescriptorsArgs
 * @property {string} smiles - SMILES string for descriptor calculation.
 */

/**
 * Calculate physicochemical properties for a molecule from its SMILES string.
 * @param {ChemblComputeMolecularDescriptorsArgs} args - The arguments for molecular descriptors calculation
 * @returns {Promise<object>} A promise that resolves to the molecular descriptors data
 */
export const computeMolecularDescriptors = async (args) => {
  if (!args || !args.smiles) {
    console.error('computeMolecularDescriptors: smiles is required');
    return Promise.reject(new Error('SMILES string is required'));
  }
  
  if (typeof args.smiles !== 'string' || args.smiles.trim() === '') {
    console.error('computeMolecularDescriptors: smiles must be a non-empty string');
    return Promise.reject(new Error('SMILES string must be a non-empty string'));
  }
  
  try {
    const result = await callMcpTool('/mcp/chembl/compute_molecular_descriptors', args);
    console.log('ChEMBL molecular descriptors computed successfully');
    return result;
  } catch (error) {
    console.error('Failed to compute molecular descriptors:', error.message);
    throw error;
  }
};

/**
 * @typedef {Object} ChemblComputeStructuralAlertsArgs
 * @property {string} smiles - SMILES string for structural alert analysis.
 */

/**
 * Identify known structural alerts within a molecule from its SMILES string.
 * @param {ChemblComputeStructuralAlertsArgs} args - The arguments for structural alerts computation
 * @returns {Promise<object>} A promise that resolves to the structural alerts data
 */
export const computeStructuralAlerts = async (args) => {
  if (!args || !args.smiles) {
    console.error('computeStructuralAlerts: smiles is required');
    return Promise.reject(new Error('SMILES string is required'));
  }
  
  if (typeof args.smiles !== 'string' || args.smiles.trim() === '') {
    console.error('computeStructuralAlerts: smiles must be a non-empty string');
    return Promise.reject(new Error('SMILES string must be a non-empty string'));
  }
  
  try {
    const result = await callMcpTool('/mcp/chembl/compute_structural_alerts', args);
    console.log('ChEMBL structural alerts computed successfully');
    return result;
  } catch (error) {
    console.error('Failed to compute structural alerts:', error.message);
    throw error;
  }
};

/**
 * @typedef {Object} ChemblStandardizeMoleculeFromSmilesArgs
 * @property {string} smiles - SMILES string of the molecule to standardize.
 */

/**
 * Standardize a molecular structure provided as a SMILES string.
 * @param {ChemblStandardizeMoleculeFromSmilesArgs} args - The arguments for molecule standardization
 * @returns {Promise<object>} A promise that resolves to the standardized molecule data
 */
export const standardizeMoleculeFromSmiles = async (args) => {
  if (!args || !args.smiles) {
    console.error('standardizeMoleculeFromSmiles: smiles is required');
    return Promise.reject(new Error('SMILES string is required'));
  }
  
  if (typeof args.smiles !== 'string' || args.smiles.trim() === '') {
    console.error('standardizeMoleculeFromSmiles: smiles must be a non-empty string');
    return Promise.reject(new Error('SMILES string must be a non-empty string'));
  }
  
  try {
    const result = await callMcpTool('/mcp/chembl/standardize_molecule_from_smiles', args);
    console.log('ChEMBL molecule standardization completed successfully');
    return result;
  } catch (error) {
    console.error('Failed to standardize molecule:', error.message);
    throw error;
  }
};

/**
 * @typedef {Object} ChemblGetParentMoleculeFromSmilesArgs
 * @property {string} smiles - SMILES string of the molecule (possibly salt/mixture).
 */

/**
 * Identify and return the parent structure for a molecule (e.g., removing counter-ions).
 * @param {ChemblGetParentMoleculeFromSmilesArgs} args - The arguments for parent molecule extraction
 * @returns {Promise<object>} A promise that resolves to the parent molecule data
 */
export const getParentMoleculeFromSmiles = async (args) => {
  if (!args || !args.smiles) {
    console.error('getParentMoleculeFromSmiles: smiles is required');
    return Promise.reject(new Error('SMILES string is required'));
  }
  
  if (typeof args.smiles !== 'string' || args.smiles.trim() === '') {
    console.error('getParentMoleculeFromSmiles: smiles must be a non-empty string');
    return Promise.reject(new Error('SMILES string must be a non-empty string'));
  }
  
  try {
    const result = await callMcpTool('/mcp/chembl/get_parent_molecule_from_smiles', args);
    console.log('ChEMBL parent molecule extraction completed successfully');
    return result;
  } catch (error) {
    console.error('Failed to get parent molecule:', error.message);
    throw error;
  }
};

// ==================== PUBCHEM TOOLS ====================

/**
 * @typedef {Object} PubchemGetCompoundByCidArgs
 * @property {string} cid - PubChem Compound ID (CID).
 */

/**
 * Retrieve compound details from PubChem using its Compound ID (CID).
 * @param {PubchemGetCompoundByCidArgs} args - The arguments for compound retrieval
 * @returns {Promise<object>} A promise that resolves to the compound data
 */
export const getCompoundByCid = async (args) => {
  if (!args || !args.cid) {
    console.error('getCompoundByCid: cid is required');
    return Promise.reject(new Error('PubChem Compound ID (CID) is required'));
  }
  
  if (typeof args.cid !== 'string' || args.cid.trim() === '') {
    console.error('getCompoundByCid: cid must be a non-empty string');
    return Promise.reject(new Error('PubChem Compound ID (CID) must be a non-empty string'));
  }
  
  try {
    const result = await callMcpTool('/mcp/pubchem/get_compound_by_cid', args);
    console.log('PubChem compound retrieval completed successfully');
    return result;
  } catch (error) {
    console.error('Failed to get PubChem compound by CID:', error.message);
    throw error;
  }
};

/**
 * @typedef {Object} PubchemGetCidsByNameArgs
 * @property {string} name - Name of the compound to search for.
 */

/**
 * Search for PubChem Compound IDs (CIDs) by compound name.
 * @param {PubchemGetCidsByNameArgs} args - The arguments for CID search by name
 * @returns {Promise<object>} A promise that resolves to the CID search results
 */
export const getCidsByName = async (args) => {
  if (!args || !args.name) {
    console.error('getCidsByName: name is required');
    return Promise.reject(new Error('Compound name is required'));
  }
  
  if (typeof args.name !== 'string' || args.name.trim() === '') {
    console.error('getCidsByName: name must be a non-empty string');
    return Promise.reject(new Error('Compound name must be a non-empty string'));
  }
  
  try {
    const result = await callMcpTool('/mcp/pubchem/get_cids_by_name', args);
    console.log('PubChem CIDs search by name completed successfully');
    return result;
  } catch (error) {
    console.error('Failed to get PubChem CIDs by name:', error.message);
    throw error;
  }
};

/**
 * @typedef {Object} PubchemGetCompoundPropertiesArgs
 * @property {string} cid - PubChem Compound ID (CID).
 * @property {string[]} properties - List of properties to retrieve.
 */

/**
 * Retrieve specific properties for a given PubChem CID.
 * @param {PubchemGetCompoundPropertiesArgs} args - The arguments for compound properties retrieval
 * @returns {Promise<object>} A promise that resolves to the compound properties data
 */
export const getCompoundProperties = async (args) => {
  if (!args || !args.cid) {
    console.error('getCompoundProperties: cid is required');
    return Promise.reject(new Error('PubChem Compound ID (CID) is required'));
  }
  
  if (!args.properties) {
    console.error('getCompoundProperties: properties is required');
    return Promise.reject(new Error('Properties list is required'));
  }
  
  if (typeof args.cid !== 'string' || args.cid.trim() === '') {
    console.error('getCompoundProperties: cid must be a non-empty string');
    return Promise.reject(new Error('PubChem Compound ID (CID) must be a non-empty string'));
  }
  
  if (!Array.isArray(args.properties) || args.properties.length === 0) {
    console.error('getCompoundProperties: properties must be a non-empty array');
    return Promise.reject(new Error('Properties must be a non-empty array'));
  }
  
  try {
    const result = await callMcpTool('/mcp/pubchem/get_compound_properties', args);
    console.log('PubChem compound properties retrieval completed successfully');
    return result;
  } catch (error) {
    console.error('Failed to get PubChem compound properties:', error.message);
    throw error;
  }
};

/**
 * @typedef {Object} PubchemGetCompoundSynonymsByCidArgs
 * @property {string} cid - PubChem Compound ID (CID).
 */

/**
 * Retrieve synonyms for a compound using its PubChem CID.
 * @param {PubchemGetCompoundSynonymsByCidArgs} args - The arguments for compound synonyms retrieval
 * @returns {Promise<object>} A promise that resolves to the compound synonyms data
 */
export const getCompoundSynonymsByCid = async (args) => {
  if (!args || !args.cid) {
    console.error('getCompoundSynonymsByCid: cid is required');
    return Promise.reject(new Error('PubChem Compound ID (CID) is required'));
  }
  
  if (typeof args.cid !== 'string' || args.cid.trim() === '') {
    console.error('getCompoundSynonymsByCid: cid must be a non-empty string');
    return Promise.reject(new Error('PubChem Compound ID (CID) must be a non-empty string'));
  }
  
  try {
    const result = await callMcpTool('/mcp/pubchem/get_compound_synonyms_by_cid', args);
    console.log('PubChem compound synonyms retrieval completed successfully');
    return result;
  } catch (error) {
    console.error('Failed to get PubChem compound synonyms by CID:', error.message);
    throw error;
  }
};

/**
 * @typedef {Object} PubchemGetCompoundImagePubchemUrlArgs
 * @property {string} cid - PubChem Compound ID (CID).
 * @property {string} [image_size='small'] - Desired image size (e.g., 'small', 'large').
 */

/**
 * Get a URL for the image of a compound from PubChem using its CID.
 * @param {PubchemGetCompoundImagePubchemUrlArgs} args - The arguments for compound image URL retrieval
 * @returns {Promise<object>} A promise that resolves to the compound image URL data
 */
export const getCompoundImagePubchemUrl = async (args) => {
  if (!args || !args.cid) {
    console.error('getCompoundImagePubchemUrl: cid is required');
    return Promise.reject(new Error('PubChem Compound ID (CID) is required'));
  }
  
  if (typeof args.cid !== 'string' || args.cid.trim() === '') {
    console.error('getCompoundImagePubchemUrl: cid must be a non-empty string');
    return Promise.reject(new Error('PubChem Compound ID (CID) must be a non-empty string'));
  }
  
  try {
    const result = await callMcpTool('/mcp/pubchem/get_compound_image_pubchem_url', args);
    console.log('PubChem compound image URL retrieval completed successfully');
    return result;
  } catch (error) {
    console.error('Failed to get PubChem compound image URL:', error.message);
    throw error;
  }
};

/**
 * @typedef {Object} PubchemGetCidsBySmilesArgs
 * @property {string} smiles - SMILES string.
 */

/**
 * Search for PubChem CIDs using a SMILES string.
 * @param {PubchemGetCidsBySmilesArgs} args - The arguments for CID search by SMILES
 * @returns {Promise<object>} A promise that resolves to the CID search results
 */
export const getCidsBySmiles = async (args) => {
  if (!args || !args.smiles) {
    console.error('getCidsBySmiles: smiles is required');
    return Promise.reject(new Error('SMILES string is required'));
  }
  
  if (typeof args.smiles !== 'string' || args.smiles.trim() === '') {
    console.error('getCidsBySmiles: smiles must be a non-empty string');
    return Promise.reject(new Error('SMILES string must be a non-empty string'));
  }
  
  try {
    const result = await callMcpTool('/mcp/pubchem/get_cids_by_smiles', args);
    console.log('PubChem CIDs search by SMILES completed successfully');
    return result;
  } catch (error) {
    console.error('Failed to get PubChem CIDs by SMILES:', error.message);
    throw error;
  }
};

/**
 * @typedef {Object} PubchemGetCidsByInchikeyArgs
 * @property {string} inchikey - InChIKey.
 */

/**
 * Search for PubChem CIDs using an InChIKey.
 * @param {PubchemGetCidsByInchikeyArgs} args - The arguments for CID search by InChIKey
 * @returns {Promise<object>} A promise that resolves to the CID search results
 */
export const getCidsByInchikey = async (args) => {
  if (!args || !args.inchikey) {
    console.error('getCidsByInchikey: inchikey is required');
    return Promise.reject(new Error('InChIKey is required'));
  }
  
  if (typeof args.inchikey !== 'string' || args.inchikey.trim() === '') {
    console.error('getCidsByInchikey: inchikey must be a non-empty string');
    return Promise.reject(new Error('InChIKey must be a non-empty string'));
  }
  
  try {
    const result = await callMcpTool('/mcp/pubchem/get_cids_by_inchikey', args);
    console.log('PubChem CIDs search by InChIKey completed successfully');
    return result;
  } catch (error) {
    console.error('Failed to get PubChem CIDs by InChIKey:', error.message);
    throw error;
  }
};

/**
 * @typedef {Object} PubchemFastIdentitySearchByCidArgs
 * @property {string} cid - PubChem Compound ID (CID) for query.
 */

/**
 * Perform a fast identity search for compounds similar to a given PubChem CID.
 * @param {PubchemFastIdentitySearchByCidArgs} args - The arguments for fast identity search
 * @returns {Promise<object>} A promise that resolves to the identity search results
 */
export const fastIdentitySearchByCid = async (args) => {
  if (!args || !args.cid) {
    console.error('fastIdentitySearchByCid: cid is required');
    return Promise.reject(new Error('PubChem Compound ID (CID) is required'));
  }
  
  if (typeof args.cid !== 'string' || args.cid.trim() === '') {
    console.error('fastIdentitySearchByCid: cid must be a non-empty string');
    return Promise.reject(new Error('PubChem Compound ID (CID) must be a non-empty string'));
  }
  
  try {
    const result = await callMcpTool('/mcp/pubchem/fast_identity_search_by_cid', args);
    console.log('PubChem fast identity search completed successfully');
    return result;
  } catch (error) {
    console.error('Failed to perform PubChem fast identity search:', error.message);
    throw error;
  }
};

/**
 * @typedef {Object} PubchemFastSubstructureSearchBySmilesArgs
 * @property {string} smiles - SMILES string for substructure query.
 */

/**
 * Perform a fast substructure search using a SMILES string.
 * @param {PubchemFastSubstructureSearchBySmilesArgs} args - The arguments for fast substructure search
 * @returns {Promise<object>} A promise that resolves to the substructure search results
 */
export const fastSubstructureSearchBySmiles = async (args) => {
  if (!args || !args.smiles) {
    console.error('fastSubstructureSearchBySmiles: smiles is required');
    return Promise.reject(new Error('SMILES string is required'));
  }
  
  if (typeof args.smiles !== 'string' || args.smiles.trim() === '') {
    console.error('fastSubstructureSearchBySmiles: smiles must be a non-empty string');
    return Promise.reject(new Error('SMILES string must be a non-empty string'));
  }
  
  try {
    const result = await callMcpTool('/mcp/pubchem/fast_substructure_search_by_smiles', args);
    console.log('PubChem fast substructure search completed successfully');
    return result;
  } catch (error) {
    console.error('Failed to perform PubChem fast substructure search:', error.message);
    throw error;
  }
};

/**
 * @typedef {Object} PubchemFastSimilarity2dSearchByCidArgs
 * @property {string} cid - PubChem Compound ID (CID) for similarity query.
 */

/**
 * Perform a fast 2D similarity search for compounds similar to a given PubChem CID.
 * @param {PubchemFastSimilarity2dSearchByCidArgs} args - The arguments for fast 2D similarity search
 * @returns {Promise<object>} A promise that resolves to the similarity search results
 */
export const fastSimilarity2dSearchByCid = async (args) => {
  if (!args || !args.cid) {
    console.error('fastSimilarity2dSearchByCid: cid is required');
    return Promise.reject(new Error('PubChem Compound ID (CID) is required'));
  }
  
  if (typeof args.cid !== 'string' || args.cid.trim() === '') {
    console.error('fastSimilarity2dSearchByCid: cid must be a non-empty string');
    return Promise.reject(new Error('PubChem Compound ID (CID) must be a non-empty string'));
  }
  
  try {
    const result = await callMcpTool('/mcp/pubchem/fast_similarity_2d_search_by_cid', args);
    console.log('PubChem fast 2D similarity search completed successfully');
    return result;
  } catch (error) {
    console.error('Failed to perform PubChem fast 2D similarity search:', error.message);
    throw error;
  }
};

/**
 * @typedef {Object} PubchemGetCidsByXrefArgs
 * @property {string} xref_type - Namespace of the cross-reference (e.g., 'RegistryID', 'RN', 'PubMedID').
 * @property {string} xref_value - Cross-reference ID.
 */

/**
 * Search for PubChem CIDs using a cross-reference ID and its namespace.
 * @param {PubchemGetCidsByXrefArgs} args - The arguments for CID search by cross-reference
 * @returns {Promise<object>} A promise that resolves to the CID search results
 */
export const getCidsByXref = async (args) => {
  if (!args || !args.xref_type) {
    console.error('getCidsByXref: xref_type is required');
    return Promise.reject(new Error('Cross-reference type is required'));
  }
  
  if (!args.xref_value) {
    console.error('getCidsByXref: xref_value is required');
    return Promise.reject(new Error('Cross-reference value is required'));
  }
  
  if (typeof args.xref_type !== 'string' || args.xref_type.trim() === '') {
    console.error('getCidsByXref: xref_type must be a non-empty string');
    return Promise.reject(new Error('Cross-reference type must be a non-empty string'));
  }
  
  if (typeof args.xref_value !== 'string' || args.xref_value.trim() === '') {
    console.error('getCidsByXref: xref_value must be a non-empty string');
    return Promise.reject(new Error('Cross-reference value must be a non-empty string'));
  }
  
  try {
    const result = await callMcpTool('/mcp/pubchem/get_cids_by_xref', args);
    console.log('PubChem CIDs search by cross-reference completed successfully');
    return result;
  } catch (error) {
    console.error('Failed to get PubChem CIDs by cross-reference:', error.message);
    throw error;
  }
};

/**
 * @typedef {Object} PubchemGetCidsByMassArgs
 * @property {string} mass_type - Type of mass to search by (e.g., 'exact', 'monoisotopic').
 * @property {number} value_or_min - Molecular mass value to search for.
 */

/**
 * Search for PubChem CIDs using molecular mass and an optional tolerance.
 * @param {PubchemGetCidsByMassArgs} args - The arguments for CID search by mass
 * @returns {Promise<object>} A promise that resolves to the CID search results
 */
export const getCidsByMass = async (args) => {
  if (!args || !args.mass_type) {
    console.error('getCidsByMass: mass_type is required');
    return Promise.reject(new Error('Mass type is required'));
  }
  
  if (args.value_or_min === undefined || args.value_or_min === null) {
    console.error('getCidsByMass: value_or_min is required');
    return Promise.reject(new Error('Mass value is required'));
  }
  
  if (typeof args.mass_type !== 'string' || args.mass_type.trim() === '') {
    console.error('getCidsByMass: mass_type must be a non-empty string');
    return Promise.reject(new Error('Mass type must be a non-empty string'));
  }
  
  if (typeof args.value_or_min !== 'number' || args.value_or_min <= 0) {
    console.error('getCidsByMass: value_or_min must be a positive number');
    return Promise.reject(new Error('Mass value must be a positive number'));
  }
  
  try {
    const result = await callMcpTool('/mcp/pubchem/get_cids_by_mass', args);
    console.log('PubChem CIDs search by mass completed successfully');
    return result;
  } catch (error) {
    console.error('Failed to get PubChem CIDs by mass:', error.message);
    throw error;
  }
};

// ==================== GENERIC MCP TOOLS EXPORTS ====================

/**
 * Generic function to call any MCP tool endpoint.
 * This is exported for flexibility in case you need to call tools not specifically implemented above.
 * @param {string} toolPath - The tool path (e.g., '/mcp/admet/get_admet_prediction')
 * @param {object} toolArgs - The arguments to pass to the tool
 * @param {string} method - HTTP method (default: 'POST')
 * @returns {Promise<object>} The response data from the tool endpoint
 */
export const callGenericMcpTool = callMcpTool;

/**
 * Ping the MCP tools server to check if it's running.
 * @returns {Promise<object>} The response data from the ping endpoint
 */
export const pingMcpToolsServer = async () => {
  try {
    const response = await mcpToolsApiClient.get('/health');
    return response.data;
  } catch (error) {
    console.error('Error pinging MCP tools server:', error.response ? error.response.data : error.message);
    throw error;
  }
};