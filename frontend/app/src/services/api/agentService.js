import axios from 'axios';
import AuthService from './AuthService.ts';

const AGENT_BASE_URL = 'http://localhost:8015'; // Base URL for your FastAPI agent

/**
 * Retrieves the auth token from localStorage via AuthService.
 */
const getAuthToken = () => {
  return AuthService.getAccessToken();
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
 * Response interceptor to handle token refresh on 401 errors.
 */
agentApiClient.interceptors.response.use(
  (response) => {
    return response;
  },
  async (error) => {
    const originalRequest = error.config;

    if (error.response?.status === 401 && !originalRequest._retry) {
      originalRequest._retry = true;

      try {
        console.log('Agent Service: Access token expired, attempting refresh...');
        const refreshSuccess = await AuthService.refreshToken();
        
        if (refreshSuccess) {
          console.log('Agent Service: Token refresh successful, retrying original request...');
          const newToken = AuthService.getAccessToken();
          
          if (newToken) {
            // Update the Authorization header for the retry
            originalRequest.headers['Authorization'] = `Bearer ${newToken}`;
            return agentApiClient.request(originalRequest);
          }
        }
      } catch (refreshError) {
        console.error('Agent Service: Token refresh failed:', refreshError);
      }

      // If refresh failed or no new token, clear auth data
      console.log('Agent Service: Token refresh failed, clearing auth data...');
      AuthService.clearAuthData();
      
      // Dispatch a custom event to notify the app about the logout
      if (typeof window !== 'undefined') {
        window.dispatchEvent(new CustomEvent('auth:logout', { 
          detail: { reason: 'token_refresh_failed', source: 'agent_service' }
        }));
      }
    }

    return Promise.reject(error);
  }
);

/**
 * Processes an SSE stream, handling events and calling appropriate callbacks.
 * @param {ReadableStreamDefaultReader} reader - The stream reader.
 * @param {object} callbacks - Callbacks for handling different stages of the SSE stream.
 * @param {(data: object) => void} [callbacks.onSessionCreated] - Called for 'session_created' event.
 * @param {(data: object) => void} callbacks.onAgentMessage - Called for 'agent_message' event.
 * @param {(data: object) => void} callbacks.onStreamEnd - Called for 'stream_end' event.
 * @param {(error: Error) => void} callbacks.onProcessingError - Called on error during stream processing.
 * @param {string} logContext - A string to identify the context in logs (e.g., "session creation", "query").
 * @returns {Promise<void>}
 */
const _processSseStream = async (reader, { onSessionCreated, onAgentMessage, onStreamEnd, onProcessingError }, logContext) => {
  const decoder = new TextDecoder();
  let buffer = '';

  try {
    while (true) {
      const { done, value } = await reader.read();
      if (done) {
        console.log(`SSE stream finished for ${logContext}.`);
        // if (buffer.trim()) { // User requested to avoid this warning.
        // console.warn(`SSE stream ended with unprocessed buffer for ${logContext}:`, buffer);
        // }
        break;
      }

      buffer += decoder.decode(value, { stream: true });
      let eolIndex;

      while ((eolIndex = buffer.indexOf('<|sep|>')) >= 0) {
        const message = buffer.substring(0, eolIndex);
        // console.log("message >>", message);
        buffer = buffer.substring(eolIndex + 7);

        let eventType = null;
        let jsonData = null;

        const lines = message.split('<|split|>');
        // console.log("lines >>", message);
        for (const line of lines) {
          // console.log("line >>", line);
          if (line.startsWith('event:')) {
            eventType = line.substring('event:'.length).trim();
          } else if (line.startsWith('data:')) {
            jsonData = line.substring('data:'.length).trim();
          }
          else {
            jsonData = line.trim();
          }
          if (eventType && jsonData) {
            try {
              const parsedData = JSON.parse(jsonData);

              if (eventType === 'session_created' && onSessionCreated) {
                onSessionCreated(parsedData);
              } else if (eventType === 'agent_message' && onAgentMessage) {
                onAgentMessage(parsedData);
              } else if (eventType === 'stream_end' && onStreamEnd) {
                onStreamEnd(parsedData);
                if (reader) await reader.cancel(); // Cancel the reader as stream has ended
                return; // Exit after stream_end
              }
            } catch (e) {
              // JSON.parse failed.
              if (eventType === 'agent_message' && onAgentMessage) {
                // If it's an agent_message, assume the jsonData is plain text (e.g., markdown).
                console.log(`Agent message is not JSON, treating as plain text (${logContext}). Data:`, jsonData);
                // Wrap the plain text in a structured object for the callback.
                // The consuming component will need to handle this structure.
                // console.log('unparsable json data >>', jsonData)
                onAgentMessage({ type: 'markdown_text', content: jsonData });
              } else {
                // For other event types (like 'session_created', 'stream_end'),
                // or if no onAgentMessage callback is provided for an agent_message,
                // this is treated as a genuine parsing error.
                const parseError = new Error(`Error parsing SSE JSON data for event '${eventType}' (${logContext}): ${e.message}. Original data: "${jsonData}"`);
                if (onProcessingError) {
                  onProcessingError(parseError);
                } else {
                  console.error(parseError.message);
                }
              }
            }
          }
        }
      }
    }
  } catch (error) {
    const streamErr = error instanceof Error ? error : new Error(String(error));
    if (onProcessingError) onProcessingError(streamErr);
    else console.error(`Error processing SSE stream for ${logContext}:`, streamErr);
    throw streamErr; // Re-throw to be caught by the caller's finally block for cleanup
  }
};

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

/**
 * @typedef {Object} SessionCreationRequestPayload
 * @property {string} query - The initial query for the session.
 */

/**
 * Creates a new user session and processes the initial query via SSE.
 * @param {SessionCreationRequestPayload} sessionRequest - The request payload containing the query.
 * @param {object} callbacks - Callbacks for handling different stages of the SSE stream.
 * @param {(data: object) => void} [callbacks.onSessionCreated] - Called when the session is created.
 * @param {(data: object) => void} [callbacks.onAgentMessage] - Called for each message from the agent.
 * @param {(data: object) => void} [callbacks.onStreamEnd] - Called when the stream ends.
 * @param {(error: Error) => void} [callbacks.onProcessingError] - Called on error during stream processing.
 * @param {(error: Error) => void} [callbacks.onSetupError] - Called on error during fetch setup.
 * @returns {Promise<void>} A promise that resolves when the stream handling is complete or rejects on error.
 */
export const createUserSession = async (
  sessionRequest,
  { onSessionCreated, onAgentMessage, onStreamEnd, onProcessingError, onSetupError }
) => {
  const token = getAuthToken();
  const headers = {
    'Content-Type': 'application/json',
  };
  if (token) {
    headers['Authorization'] = `Bearer ${token}`;
  }

  let reader; // Declare reader here to access in finally block

  try {
    const response = await fetch(`${AGENT_BASE_URL}/sessions`, {
      method: 'POST',
      headers: headers,
      body: JSON.stringify(sessionRequest),
    });

    if (!response.ok) {
      const errorData = await response.json().catch(() => ({ detail: response.statusText }));
      const error = new Error(`HTTP error! status: ${response.status} - ${errorData.detail || response.statusText}`);
      if (onSetupError) onSetupError(error);
      else console.error('Error creating user session (initial response):', errorData);
      throw error;
    }

    if (!response.body) {
      const error = new Error('Response body is null');
      if (onSetupError) onSetupError(error);
      else console.error(error.message);
      throw error;
    }

    reader = response.body.getReader();
    console.log('Starting to process SSE stream for session creation...');

    // This outer try/catch handles errors in stream processing
    try {
      await _processSseStream(reader, { onSessionCreated, onAgentMessage, onStreamEnd, onProcessingError }, 'session creation');
    } catch (error) {
      // Error already logged by _processSseStream, re-throw if needed or handle
      // console.error('Error during _processSseStream for session creation:', error);
      // No need to call onProcessingError here as _processSseStream already does.
      throw error; // Re-throw to ensure finally block cleans up
    } finally {
      if (reader) {
        // Ensure the lock is released if the stream wasn't fully consumed or cancelled.
        // reader.cancel() should ideally handle this, but as a safeguard:
        try {
          if (reader.locked) { // Check if locked before trying to release.
            reader.releaseLock();
          }
        } catch (releaseLockError) {
          // console.warn("Error releasing reader lock:", releaseLockError);
          // This can happen if cancel was already called and released it.
        }
      }
    }

  } catch (error) {
    const setupErr = error instanceof Error ? error : new Error(String(error));
    if (onSetupError) onSetupError(setupErr);
    else console.error('Error creating user session (fetch setup):', setupErr.message);
    throw setupErr;
  }
};

/**
 * Lists sessions for the current user. If sessionId is provided, fetches a specific session.
 * @returns {Promise<Object>} The response data from the endpoint.
 */
export const listUserSessions = async (sessionId) => {
  try {
    const url = `/sessions`;
    const token = getAuthToken();
    if (!token) {
      console.error('No token found for listUserSessions');
      return Promise.reject(new Error('Authentication token not found'));
    }
    const response = await agentApiClient.get(url,
      {
        method: 'GET',
        headers: {
          'Authorization': `Bearer ${token}`,
          'Content-Type': 'application/json'
        },
      }
    );
    return response.data;
  } catch (error) {
    console.error('Error listing user sessions:', error.response ? error.response.data : error.message);
    if (error.response && error.response.status === 401) {
      console.error('Authentication failed. Token might be invalid or expired.');
    }
    throw error;
  }
};

/**
 * Fetches details for a specific session.
 * @param {string} sessionId The ID of the session to fetch.
 * @returns {Promise<object>} A promise that resolves to the session details.
 */
export const getSessionDetails = async (sessionId) => {
  if (!sessionId) {
    console.error('getSessionDetails: sessionId is required');
    return Promise.reject(new Error('Session ID is required'));
  }
  const token = getAuthToken(); 
  if (!token) {
    console.error('No token found for getSessionDetails');
    return Promise.reject(new Error('Authentication token not found'));
  }

  try {
    const response = await fetch(`${AGENT_BASE_URL}/sessions/${sessionId}`, {
      method: 'GET',
      headers: {
        'Authorization': `Bearer ${token}`,
        'Content-Type': 'application/json'
      },
    });

    if (!response.ok) {
      const errorData = await response.text();
      console.error('Error fetching session details:', response.status, errorData);
      throw new Error(`HTTP error ${response.status}: ${errorData}`);
    }
    return await response.json();
  } catch (error) {
    console.error('Failed to fetch session details:', error);
    throw error;
  }
};

/**
 * Sends a follow-up query to an existing session.
 * This will be an SSE endpoint.
 * @param {string} sessionId The ID of the session.
 * @param {string} queryText The query text to send.
 * @param {object} callbacks Callbacks for SSE events.
 * @param {(data: any) => void} callbacks.onAgentMessage Callback for agent messages.
 * @param {(data: any) => void} callbacks.onStreamEnd Callback for stream end.
 * @param {(error: Error) => void} callbacks.onProcessingError Callback for errors during SSE processing.
 * @param {(error: Error) => void} callbacks.onSetupError Callback for errors during SSE setup (e.g., network issues).
 */
export const sendQueryToSession = async (sessionId, queryText, { onAgentMessage, onStreamEnd, onProcessingError, onSetupError }) => {
  if (!sessionId || !queryText) {
    console.error('sendQueryToSession: sessionId and queryText are required');
    if (onSetupError) onSetupError(new Error('Session ID and query text are required'));
    return;
  }
  
  const token = getAuthToken();
  if (!token) {
    console.error('No token found for sendQueryToSession');
    if (onSetupError) onSetupError(new Error('Authentication token not found'));
    return;
  }

  let reader;
  try {
    const response = await fetch(`${AGENT_BASE_URL}/sessions/${sessionId}/query`, {
      method: 'POST',
      headers: {
        'Authorization': `Bearer ${token}`,
        'Content-Type': 'application/json',
        'Accept': 'text/event-stream'
      },
      body: JSON.stringify({ query: queryText }), // Changed to match API's SessionCreationRequest format
    });

    if (!response.ok) {
      const errorText = await response.text();
      const setupError = new Error(`HTTP error ${response.status}: ${errorText}`);
      console.error('Error sending query to session:', setupError.message);
      if (onSetupError) onSetupError(setupError);
      return;
    }

    if (!response.body) {
      const setupError = new Error('Response body is null');
      if (onSetupError) onSetupError(setupError);
      return;
    }
    reader = response.body.getReader();
    console.log('Starting to process SSE stream for query...');

    try {
      await _processSseStream(reader, { onAgentMessage, onStreamEnd, onProcessingError }, 'query');
    } catch (error) {
      // Error already logged by _processSseStream.
      // No need to call onProcessingError here as _processSseStream already does.
      // The main catch block for sendQueryToSession will handle onSetupError if this was a setup issue,
      // or this re-thrown error will be caught by it.
      // console.error('Error during _processSseStream for query:', error);
      throw error; // Re-throw to ensure finally block cleans up
    } finally {
      if (reader) {
        // Ensure the lock is released if the stream wasn't fully consumed or cancelled.
        // reader.cancel() should ideally handle this, but as a safeguard:
        try {
          if (reader.locked) { // Check if locked before trying to release.
            reader.releaseLock();
          }
        } catch (releaseLockError) {
          // console.warn("Error releasing reader lock (query):", releaseLockError);
        }
      }
    }

  } catch (error) {
    console.error('Failed to send query to session:', error);
    if (onSetupError) onSetupError(error);
  }
};

/**
 * Retrieves a complete tool response by its ID.
 * @param {string} toolId - The unique identifier of the tool response to retrieve.
 * @returns {Promise<object>} A promise that resolves to the complete tool response data.
 */
export const getToolResponse = async (toolId) => {
  if (!toolId) {
    console.error('getToolResponse: toolId is required');
    return Promise.reject(new Error('Tool ID is required'));
  }

  const token = getAuthToken();
  if (!token) {
    console.error('No token found for getToolResponse');
    return Promise.reject(new Error('Authentication token not found'));
  }
  
  try {
    // Using fetch instead of axios to match the pattern in createUserSession
    const response = await fetch(`${AGENT_BASE_URL}/tool-responses/${toolId}`, {
      method: 'GET',
      headers: {
        'Authorization': `Bearer ${token}`,
        'Content-Type': 'application/json'
      },
    });

    if (!response.ok) {
      const errorData = await response.text();
      console.error('Error fetching tool response:', response.status, errorData);
      throw new Error(`HTTP error ${response.status}: ${errorData}`);
    }
    
    return await response.json();
  } catch (error) {
    console.error('Error retrieving tool response:', error.message);
    throw error;
  }
};

/**
 * Retrieves all tool responses for a specific session.
 * @param {string} sessionId - The unique identifier of the session.
 * @returns {Promise<object>} A promise that resolves to an object containing all tool responses for the session.
 */
export const getSessionToolResponses = async (sessionId) => {
  if (!sessionId) {
    console.error('getSessionToolResponses: sessionId is required');
    return Promise.reject(new Error('Session ID is required'));
  }
  
  try {
    const response = await agentApiClient.get(`/sessions/${sessionId}/tool-responses`);
    return response.data;
  } catch (error) {
    console.error('Error retrieving session tool responses:', error.response ? error.response.data : error.message);
    if (error.response && error.response.status === 404) {
      console.error(`Session with ID ${sessionId} not found or does not belong to the current user`);
    }
    throw error;
  }
};

/**
 * Retrieves the event state of a specific session using its ID.
 * @param {string} sessionId - The unique identifier of the session.
 * @returns {Promise<object>} A promise that resolves to the session's event state.
 */
export const getSessionEventState = async (sessionId) => {
  if (!sessionId) {
    console.error('getSessionEventState: sessionId is required');
    return Promise.reject(new Error('Session ID is required'));
  }

  const token = getAuthToken();
  if (!token) {
    console.error('No token found for getSessionEventState');
    return Promise.reject(new Error('Authentication token not found'));
  }

  try {
    const response = await fetch(`${AGENT_BASE_URL}/sessions/${sessionId}`, {
      method: 'GET',
      headers: {
        'Authorization': `Bearer ${token}`,
        'Content-Type': 'application/json'
      },
    });

    if (!response.ok) {
      const errorData = await response.text();
      console.error('Error fetching session event state:', response.status, errorData);
      throw new Error(`HTTP error ${response.status}: ${errorData}`);
    }
    
    const data = await response.json();
    if (data) {
      // console.log('Retrieved session event state:', data.state);
      // return data.state;
    // } else {
      // console.warn('Session retrieved but no state field found');
      return data; // Return the full response if state isn't specifically available
    }
  } catch (error) {
    console.error('Error retrieving session event state:', error.message);
    throw error;
  }
};

/**
 * Deletes a specific session using its ID.
 * @param {string} sessionId - The unique identifier of the session to delete.
 * @returns {Promise<object>} A promise that resolves to the response data from the deletion.
 */
export const deleteSession = async (sessionId) => {
  if (!sessionId) {
    console.error('deleteSession: sessionId is required');
    return Promise.reject(new Error('Session ID is required'));
  }

  const token = getAuthToken();
  if (!token) {
    console.error('No token found for deleteSession');
    return Promise.reject(new Error('Authentication token not found'));
  }

  try {
    const response = await fetch(`${AGENT_BASE_URL}/sessions/${sessionId}`, {
      method: 'DELETE',
      headers: {
        'Authorization': `Bearer ${token}`,
        'Content-Type': 'application/json'
      },
    });

    if (!response.ok) {
      const errorData = await response.text();
      console.error('Error deleting session:', response.status, errorData);
      throw new Error(`HTTP error ${response.status}: ${errorData}`);
    }
    
    return await response.json();
  } catch (error) {
    console.error('Error deleting session:', error.message);
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
