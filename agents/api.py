from fastapi import FastAPI, Depends, Request, HTTPException, Query
from fastapi.responses import StreamingResponse, JSONResponse # Added for SSE and JSON responses
from fastapi.middleware.cors import CORSMiddleware
from contextlib import asynccontextmanager
import json # Added for SSE data formatting

from google.adk.sessions import DatabaseSessionService
from google.adk.runners import Runner 

from auth_utils import verify_token
from master_agent.agent import root_agent as master_agent 
from run_utils import call_agent_async 
from custom_utils.mongodb_handler import MongoDBHandler  # Import MongoDB handler

import uvicorn
import os
from dotenv import load_dotenv
from pydantic import BaseModel, Field

AGENT_APP_NAME = "Master Agent"

# Lifespan manager for application setup/teardown
@asynccontextmanager
async def lifespan(app: FastAPI):
    # Load resources or setup agents here
    print("Starting up agent services...")
    # Initialize Database Session Service
    db_url = os.getenv("DATABASE_URL", "sqlite:///./agent_api_data.db")
    app.state.session_service = DatabaseSessionService(db_url=db_url)
    print(f"Database session service initialized with URL: {db_url}")
    
    # Initialize MongoDB connection
    try:
        app.state.mongodb_handler = MongoDBHandler()
        print("MongoDB handler initialized")
    except Exception as e:
        print(f"Error initializing MongoDB handler: {str(e)}")
    
    yield
    # Clean up resources or teardown agents here
    print("Shutting down agent services...")


app = FastAPI(lifespan=lifespan)

load_dotenv()

# CORS configuration
# Adjust origins based on where your frontend/clients are hosted
origins = [
    "http://localhost",
    "http://localhost:3000",
    "http://localhost:5173",
    "http://localhost:8080",
    "http://127.0.0.1",
    "http://127.0.0.1:3000",
    "http://127.0.0.1:5173",
    "http://127.0.0.1:8080",
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.get("/ping")
async def ping():
    return {"message": "pong"}


@app.get("/test")
async def test_endpoint(current_user: dict = Depends(verify_token)):
    return {"message": "Test endpoint reached successfully", "user_details": current_user}


class SessionCreationRequest(BaseModel):
    query: str


async def _event_generator_for_new_session(
    session_service: DatabaseSessionService,
    current_user_payload: dict,
    session_creation_req: SessionCreationRequest,
    agent_module, 
    app_name_str: str
):
    user_id = current_user_payload.get("id")
    user_name = current_user_payload.get("username", current_user_payload.get("name", "Unknown User"))

    # ADK's create_session is synchronous, remove await
    new_session = await session_service.create_session(
        app_name=app_name_str,
        user_id=str(user_id),
        state={"user_name": user_name} # Add initial state
    )

    session_created_data = {
        "session_id": new_session.id,
        "message": f"New session created for user '{user_name}' (ID: {user_id}). Processing initial query..."
    }
    # json.dumps for session_created_data will be a single line by default.
    yield f"event: session_created<|split|>data: {json.dumps(session_created_data)}<|sep|>"

    runner = Runner(
        agent=agent_module,
        app_name=app_name_str,
        session_service=session_service,
    )

    async for content_part in call_agent_async(
        runner, str(user_id), new_session.id, session_creation_req.query
    ):
        if content_part == "<|done|>":
            print("Query response complete, ending stream.")
            stream_end_data = {"message": "Initial query processing complete."}
            print("Sent Data >>",f"event: stream_end<|split|>data: {json.dumps(stream_end_data)}<|sep|>")
            yield f"event: stream_end<|split|>data: {json.dumps(stream_end_data)}<|sep|>"
            break
        
        # Assuming content_part is now a dictionary from run_utils.py
        if isinstance(content_part, dict):
            try:
                # # Check if this is a tool response with a tool_id
                # if "content" in content_part and isinstance(content_part["content"], str) and "tool_id:" in content_part["content"]:
                #     # Extract tool_id if present
                #     import re
                #     tool_id_match = re.search(r'tool_id: ([a-zA-Z0-9-]+)', content_part["content"])
                #     if tool_id_match:
                #         tool_id = tool_id_match.group(1)
                #         # Add tool_id to the response payload
                #         content_part["tool_id"] = tool_id
                
                # json.dumps without indent produces a single line
                json_payload = json.dumps(content_part) 
                print("Sent Data >>", f"event: agent_message<|split|>data: {json_payload}<|sep|>")
                yield f"event: agent_message<|split|>data: {json_payload}<|sep|>"
            except TypeError as e:
                error_payload = json.dumps({"type": "error", "content": f"Failed to serialize agent message: {str(e)}"})
                print("Sent Data >>", f"event: agent_message<|split|>data: {error_payload}<|sep|>")
                yield f"event: agent_message<|split|>data: {error_payload}<|sep|>"
        else:
            # Fallback for any unexpected string content (should ideally not happen if run_utils is consistent)
            unknown_payload = json.dumps({"type": "unknown_string_content", "content": str(content_part)})
            print("Sent Data >>", f"event: agent_message<|split|>data: {unknown_payload}<|sep|>")
            yield f"event: agent_message<|split|>data: {unknown_payload}<|sep|>"


@app.post("/sessions", status_code=200) # Changed status_code to 200 for StreamingResponse
async def create_user_session(
    request: Request,
    session_request: SessionCreationRequest, 
    current_user: dict = Depends(verify_token)
):
    return StreamingResponse(
        _event_generator_for_new_session(
            session_service=request.app.state.session_service,
            current_user_payload=current_user,
            session_creation_req=session_request,
            agent_module=master_agent, 
            app_name_str=AGENT_APP_NAME
        ),
        media_type="text/event-stream"
    )

@app.get("/sessions")
async def list_user_sessions(
    request: Request,
    current_user: dict = Depends(verify_token)
):
    session_service: DatabaseSessionService = request.app.state.session_service
    user_id = current_user.get("id")

    if not user_id:
        raise HTTPException(status_code=400, detail="User ID not found in token")

    existing_sessions = await session_service.list_sessions(
        app_name=AGENT_APP_NAME,
        user_id=str(user_id),
    )
    return existing_sessions


@app.get("/sessions/{session_id}", status_code=200)
async def get_session_details(
    session_id: str,
    request: Request,
    current_user: dict = Depends(verify_token)
):
    session_service: DatabaseSessionService = request.app.state.session_service
    user_id = current_user.get("id")

    if not user_id:
        raise HTTPException(status_code=400, detail="User ID not found in token")

    try:
        retrieved_session = await session_service.get_session(
            session_id=session_id,
            app_name=AGENT_APP_NAME,
            user_id=str(user_id)
        )
        return retrieved_session
    except Exception as e:
        print(f"Error retrieving session {session_id} for user {user_id}: {e}")
        raise HTTPException(status_code=500, detail="Internal server error while retrieving session.")


async def _event_generator_for_session_query(
    session_service: DatabaseSessionService,
    current_user_payload: dict,
    session_id: str,
    query: str,
    agent_module,
    app_name_str: str
):
    user_id = current_user_payload.get("id")
    
    # Verify session exists and belongs to the user (optional, but good practice)
    try:
        # ADK's get_session is asynchronous
        await session_service.get_session(
            session_id=session_id,
            app_name=app_name_str,
            user_id=str(user_id)
        )
    except Exception as e:
        error_message = f"Failed to retrieve session {session_id} for user {user_id}: {str(e)}"
        print(error_message)
        error_payload = json.dumps({"type": "error", "content": error_message})
        yield f"event: agent_message\\\\ndata: {error_payload}\\\\n\\\\n"
        stream_end_data = {"message": "Query processing failed due to session error."}
        yield f"event: stream_end\\\\ndata: {json.dumps(stream_end_data)}\\\\n\\\\n"
        return

    runner = Runner(
        agent=agent_module,
        app_name=app_name_str,
        session_service=session_service,
    )

    async for content_part in call_agent_async(
        runner, str(user_id), session_id, query
    ):
        if content_part == "<|done|>":
            print(f"Query response complete for session {session_id}, ending stream.")
            stream_end_data = {"message": "Query processing complete."}
            yield f"event: stream_end\\\\ndata: {json.dumps(stream_end_data)}\\\\n\\\\n"
            break
        
        if isinstance(content_part, dict):
            try:
                # Check if this is a tool response with a tool_id
                if "content" in content_part and isinstance(content_part["content"], str) and "tool_id:" in content_part["content"]:
                    # Extract tool_id if present
                    import re
                    tool_id_match = re.search(r'tool_id: ([a-zA-Z0-9-]+)', content_part["content"])
                    if tool_id_match:
                        tool_id = tool_id_match.group(1)
                        # Add tool_id to the response payload
                        content_part["tool_id"] = tool_id
                
                json_payload = json.dumps(content_part)
                yield f"event: agent_message\\\\ndata: {json_payload}\\\\n\\\\n"
            except TypeError as e:
                error_payload = json.dumps({"type": "error", "content": f"Failed to serialize agent message: {str(e)}"})
                yield f"event: agent_message\\\\ndata: {error_payload}\\\\n\\\\n"
        else:
            unknown_payload = json.dumps({"type": "unknown_string_content", "content": str(content_part)})
            yield f"event: agent_message\\\\ndata: {unknown_payload}\\\\n\\\\n"


@app.post("/sessions/{session_id}/query", status_code=200)
async def query_session(
    session_id: str,
    request: Request,
    session_query: SessionCreationRequest, 
    current_user: dict = Depends(verify_token)
):
    # ADK's get_session is asynchronous, ensure it's awaited if called directly here for validation
    # However, the validation is now part of the _event_generator_for_session_query
    return StreamingResponse(
        _event_generator_for_session_query(
            session_service=request.app.state.session_service,
            current_user_payload=current_user,
            session_id=session_id,
            query=session_query.query,
            agent_module=master_agent,
            app_name_str=AGENT_APP_NAME
        ),
        media_type="text/event-stream"
    )


# Define a Pydantic model for tool response request
class ToolResponseRequest(BaseModel):
    tool_id: str = Field(..., description="Unique identifier for the tool response")


@app.get("/tool-responses/{tool_id}")
async def get_tool_response(
    tool_id: str,
    request: Request,
    current_user: dict = Depends(verify_token)
):
    """
    Retrieve a complete tool response by its ID.
    
    This endpoint allows retrieving the full, non-truncated tool response that was
    stored in MongoDB during tool execution.
    """
    mongodb_handler = request.app.state.mongodb_handler
    if not mongodb_handler:
        raise HTTPException(status_code=500, detail="MongoDB handler not available")
    
    # Retrieve the tool response from MongoDB
    tool_response = mongodb_handler.get_tool_response(tool_id)
    
    if not tool_response:
        raise HTTPException(status_code=404, detail=f"Tool response with ID {tool_id} not found")
    
    # Check if the user is authorized to access this tool response
    # response_user_id = tool_response.get("user_id", "unknown")
    # current_user_id = current_user.get("id")
    
    # if str(response_user_id) != str(current_user_id):
    #     raise HTTPException(status_code=403, detail="You are not authorized to access this tool response")
    
    # Return the complete tool response
    return JSONResponse(content={
        "tool_id": tool_id,
        "tool_name": tool_response.get("tool_name"),
        "session_id": tool_response.get("session_id"),
        "created_at": tool_response.get("created_at").isoformat() if tool_response.get("created_at") else None,
        "response_data": tool_response.get("response_data")
    })


@app.get("/sessions/{session_id}/tool-responses")
async def get_session_tool_responses(
    session_id: str,
    request: Request,
    current_user: dict = Depends(verify_token)
):
    """
    Retrieve all tool responses for a specific session.
    
    This endpoint allows retrieving all tool responses that were stored in MongoDB
    for a specific session.
    """
    # Verify the session exists and belongs to the user
    session_service: DatabaseSessionService = request.app.state.session_service
    user_id = current_user.get("id")
    
    try:
        await session_service.get_session(
            session_id=session_id,
            app_name=AGENT_APP_NAME,
            user_id=str(user_id)
        )
    except Exception as e:
        raise HTTPException(status_code=404, detail=f"Session {session_id} not found or does not belong to the current user")
    
    # Retrieve tool responses from MongoDB
    mongodb_handler = request.app.state.mongodb_handler
    if not mongodb_handler:
        raise HTTPException(status_code=500, detail="MongoDB handler not available")
    
    tool_responses = mongodb_handler.get_tool_responses_by_session(session_id)
    
    # Format the responses for JSON serialization
    formatted_responses = []
    for response in tool_responses:
        formatted_responses.append({
            "tool_id": response.get("tool_id"),
            "tool_name": response.get("tool_name"),
            "created_at": response.get("created_at").isoformat() if response.get("created_at") else None,
            "response_data": response.get("response_data")
        })
    
    return JSONResponse(content={"session_id": session_id, "tool_responses": formatted_responses})


if __name__ == "__main__":
    # load_dotenv() is called at module level, or uvicorn can handle .env with --env-file
    port = int(os.getenv("PORT", "8015")) # Ensure PORT is a string for getenv
    uvicorn.run(app, host="0.0.0.0", port=port)