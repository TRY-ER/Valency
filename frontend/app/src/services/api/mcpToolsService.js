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

// ==================== FUTURE TOOL CATEGORIES ====================
// Add other MCP tool categories here as they become available
// For example:
// - Molecular property prediction tools
// - Protein analysis tools  
// - Drug discovery tools
// - Chemical reaction prediction tools
// etc.

// Example usage:
/*
import { 
  getAdmetPrediction, 
  getAlphafoldPrediction,
  getUniprotSummary,
  getAlphafoldAnnotations,
  getBricsCandidates,
  callGenericMcpTool 
} from './services/api/mcpToolsService';

// Using the specific ADMET function
const handleAdmetPrediction = async () => {
  try {
    const result = await getAdmetPrediction({ 
      smiles: 'CCO' // Ethanol SMILES
    });
    console.log('ADMET Prediction:', result);
  } catch (error) {
    console.error('Error:', error.message);
  }
};

// Using the AlphaFold prediction function
const handleAlphafoldPrediction = async () => {
  try {
    const result = await getAlphafoldPrediction({ 
      uniprot_accession: 'P04637' // Example UniProt ID for p53
    });
    console.log('AlphaFold Prediction:', result);
  } catch (error) {
    console.error('Error:', error.message);
  }
};

// Using the UniProt summary function
const handleUniprotSummary = async () => {
  try {
    const result = await getUniprotSummary({ 
      qualifier: 'P04637', // UniProtKB accession number, entry name, or CRC64 checksum
      start: 1,
      end: 100
    });
    console.log('UniProt Summary:', result);
  } catch (error) {
    console.error('Error:', error.message);
  }
};

// Using the AlphaFold annotations function
const handleAlphafoldAnnotations = async () => {
  try {
    const result = await getAlphafoldAnnotations({ 
      qualifier: 'P04637',
      annotation_type: 'MUTAGEN'
    });
    console.log('AlphaFold Annotations:', result);
  } catch (error) {
    console.error('Error:', error.message);
  }
};

// Using the BRICS candidates function
const handleBricsCandidates = async () => {
  try {
    const result = await getBricsCandidates({ 
      smiles_list: ['CCO', 'C1=CC=CC=C1', 'CC(=O)O'], // Example SMILES strings
      is_polymer: false
    });
    console.log('BRICS Candidates:', result);
  } catch (error) {
    console.error('Error:', error.message);
  }
};

// Using the generic function for flexibility
const handleGenericToolCall = async () => {
  try {
    const result = await callGenericMcpTool(
      '/mcp/admet/get_admet_prediction',
      { smiles: 'CCO' }
    );
    console.log('Tool Result:', result);
  } catch (error) {
    console.error('Error:', error.message);
  }
};
*/
