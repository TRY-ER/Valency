import sys
import os

# Add the project root to sys.path to allow importing 'engine'
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

from env_utils import set_env_vars
set_env_vars("./.env")
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from app.routers import (
    base,
    validator,
    visualize,
    informer,
    generator,
    discriminator,
    chat,
)
from app.routers.mcp_adapter_router import all_mcp_routers # Import the list of routers
import uvicorn

app = FastAPI(title="Chemistry API", description="API for Chemistry", version="0.1.0")

# allow origins
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Include routers
app.include_router(base.router, tags=["base"]) # Include the base router
app.include_router(validator.router, prefix="/validate", tags=["validator"]) # Include the validator router
app.include_router(visualize.router, prefix="/visualize", tags=["visualize"]) # Include the validator router
app.include_router(informer.router, prefix="/inform", tags=["inform"]) # Include the validator router
app.include_router(generator.router, prefix="/generate", tags=["generate"]) # Include the validator router
app.include_router(discriminator.router, prefix="/discriminator", tags=["discriminate"]) # Include the validator router
app.include_router(chat.router, prefix="/chat", tags=["chat"]) # Include the validator router

# Include all MCP routers
for router_item in all_mcp_routers:
    app.include_router(router_item)

if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=8000) # Run the app