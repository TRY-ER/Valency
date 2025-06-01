import axios from 'axios';

const AGENT_BASE_URL = 'http://localhost:8015'; // Base URL for your FastAPI agent

/**
 * Retrieves the auth token from localStorage.
 * In a real app, this might come from a more robust auth service or context.
 */
const getAuthToken = () => {
  return localStorage.getItem('access_token');
};

/**
 * Creates an axios instance with common configurations for the agent service.
 */
const agentApiClient = axios.create({
  baseURL: AGENT_BASE_URL,
});

/**
 * Interceptor to add the Authorization header to requests if a token exists.
 */
agentApiClient.interceptors.request.use(
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
 * Calls the /ping endpoint.
 * @returns {Promise<Object>} The response data from the endpoint.
 */
export const pingAgent = async () => {
  try {
    const response = await agentApiClient.get('/ping');
    return response.data;
  } catch (error) {
    console.error('Error calling ping endpoint:', error.response ? error.response.data : error.message);
    throw error;
  }
};

/**
 * Calls the /test endpoint, which requires authentication.
 * @returns {Promise<Object>} The response data from the endpoint.
 */
export const testAgentEndpoint = async () => {
  try {
    const response = await agentApiClient.get('/test');
    return response.data;
  } catch (error) {
    console.error('Error calling test endpoint:', error.response ? error.response.data : error.message);
    // Handle specific auth errors if needed, e.g., redirect to login
    if (error.response && error.response.status === 401) {
      console.error('Authentication failed. Token might be invalid or expired.');
      // Optionally, trigger a logout or token refresh mechanism here
    }
    throw error;
  }
};

// Example of how you might use these functions in a component:
/*
import React, { useEffect, useState } from 'react';
import { pingAgent, testAgentEndpoint } from './services/api/agentService'; // Adjust path as needed

function AgentTester() {
  const [pingResponse, setPingResponse] = useState(null);
  const [testResponse, setTestResponse] = useState(null);
  const [error, setError] = useState('');

  useEffect(() => {
    pingAgent()
      .then(data => setPingResponse(data))
      .catch(err => setError('Ping failed: ' + err.message));

    testAgentEndpoint()
      .then(data => setTestResponse(data))
      .catch(err => setError('Test endpoint failed: ' + err.message));
  }, []);

  return (
    <div>
      <h2>Agent API Test</h2>
      {error && <p style={{ color: 'red' }}>Error: {error}</p>}
      <div>
        <h3>Ping Response:</h3>
        <pre>{JSON.stringify(pingResponse, null, 2)}</pre>
      </div>
      <div>
        <h3>Test Endpoint Response (Authenticated):</h3>
        <pre>{JSON.stringify(testResponse, null, 2)}</pre>
      </div>
    </div>
  );
}

export default AgentTester;
*/
