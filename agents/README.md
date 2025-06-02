# Valency Agents API

This is the API service for Valency Agents, providing agent session management and integration with MongoDB for storing complete tool responses.

## Setup

1. Install dependencies:
   ```bash
   pip install -e .
   ```

2. Set up environment variables:
   ```bash
   # Create a .env file
   touch .env

   # Add required variables
   echo "DATABASE_URL=sqlite:///./agent_api_data.db" >> .env
   echo "MONGODB_URI=mongodb://localhost:27017/" >> .env
   echo "MONGODB_DB_NAME=valency_agents" >> .env
   ```

3. Run MongoDB:
   Make sure you have MongoDB running locally or update the `MONGODB_URI` variable with your MongoDB connection string.

   ```bash
   # Start MongoDB using Docker (if not installed locally)
   docker run -d -p 27017:27017 --name mongodb mongo:latest
   ```

4. Start the API server:
   ```bash
   python api.py
   ```

## Tool Response Storage and Retrieval

The API includes functionality to store complete tool responses in MongoDB and retrieve them when needed:

1. When a tool response exceeds the maximum length (defined in `tool_output_overflow_callback.py`), the full response is stored in MongoDB and a truncated version is sent to the client with a tool_id reference.

2. The client can then use the tool_id to retrieve the complete response via the `/tool-responses/{tool_id}` endpoint.

3. Session tool responses can be listed via the `/sessions/{session_id}/tool-responses` endpoint.

## API Endpoints

### Session Management
- `GET /ping` - Health check endpoint
- `POST /sessions` - Create a new session
- `GET /sessions` - List user sessions
- `GET /sessions/{session_id}` - Get session details
- `POST /sessions/{session_id}/query` - Query an existing session

### Tool Response Management
- `GET /tool-responses/{tool_id}` - Get a complete tool response by ID
- `GET /sessions/{session_id}/tool-responses` - List all tool responses for a session

## Environment Variables

- `DATABASE_URL` - URL for the SQLite database (default: `sqlite:///./agent_api_data.db`)
- `MONGODB_URI` - MongoDB connection string (default: `mongodb://localhost:27017/`)
- `MONGODB_DB_NAME` - MongoDB database name (default: `valency_agents`)
- `PORT` - Port for the API server (default: `8015`)