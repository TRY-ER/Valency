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
import { getAdmetPrediction, callGenericMcpTool } from './services/api/mcpToolsService';

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
