from fastapi import FastAPI, Depends
from fastapi.middleware.cors import CORSMiddleware  # Added CORS
from contextlib import asynccontextmanager  # Added for lifespan
from .auth_utils import verify_token  # Import the new dependency
import uvicorn
import os
from dotenv import load_dotenv

# Lifespan manager for application setup/teardown
@asynccontextmanager
async def lifespan(app: FastAPI):
    # Load resources or setup agents here
    print("Starting up agent services...")
    # Example: Initialize agent managers, load configurations, etc.
    # app.state.my_agent_manager = await initialize_agent_manager()
    yield
    # Clean up resources or teardown agents here
    print("Shutting down agent services...")
    # Example: Gracefully shutdown agent connections
    # if hasattr(app.state, 'my_agent_manager'):
    #     await app.state.my_agent_manager.shutdown()


app = FastAPI(lifespan=lifespan)  # Added lifespan

load_dotenv()  # Load .env file at the module level if not already done by uvicorn --env-file

# CORS configuration
# Adjust origins based on where your frontend/clients are hosted
origins = [
    "http://localhost",         # Allow local development
    "http://localhost:3000",    # Common React/Next.js dev port
    "http://localhost:5173",    # Common Vite/Svelte dev port
    "http://localhost:8080",    # Another common dev port
    "http://127.0.0.1",
    "http://127.0.0.1:3000",
    "http://127.0.0.1:5173",
    "http://127.0.0.1:8080",
    # Add your production frontend URL here, e.g., "https://your-agent-frontend.com"
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins, # A list of specific origins allowed
    allow_credentials=True, # Support cookies
    allow_methods=["*"],    # Allow all standard methods (GET, POST, etc.)
    allow_headers=["*"],    # Allow all headers
)


@app.get("/ping")
async def ping():
    return {"message": "pong"}


@app.get("/test")
async def test_endpoint(current_user: dict = Depends(verify_token)):
    return {"message": "Test endpoint reached successfully", "user_details": current_user}


if __name__ == "__main__":
    # load_dotenv() is called at module level, or uvicorn can handle .env with --env-file
    port = int(os.getenv("PORT", "8015")) # Ensure PORT is a string for getenv
    uvicorn.run(app, host="0.0.0.0", port=port)