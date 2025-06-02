from fastapi import FastAPI, Depends, Request, HTTPException, Query
from fastapi.responses import StreamingResponse # Added for SSE
from fastapi.middleware.cors import CORSMiddleware
from contextlib import asynccontextmanager
import json # Added for SSE data formatting

from google.adk.sessions import DatabaseSessionService
from google.adk.runners import Runner 

from auth_utils import verify_token
from master_agent.agent import root_agent as master_agent 
from run_utils import call_agent_async 

import uvicorn
import os
from dotenv import load_dotenv
from pydantic import BaseModel

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
    yield f"event: session_created\\ndata: {json.dumps(session_created_data)}\\n\\n"

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


if __name__ == "__main__":
    # load_dotenv() is called at module level, or uvicorn can handle .env with --env-file
    port = int(os.getenv("PORT", "8015")) # Ensure PORT is a string for getenv
    uvicorn.run(app, host="0.0.0.0", port=port)