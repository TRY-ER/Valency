from custom_utils.repo.declare import ADMET_INSTRUCTIONS_MD_PATH # Changed import
from custom_utils import after_tool_output_limit_callback
import os
from google.adk.tools.mcp_tool.mcp_toolset import (
    MCPToolset,
    SseServerParams
)
from google.adk.agents import Agent # google.adk.agents.Agent instead of LoopAgent
from custom_utils.file_reader import read_markdown_file
from dotenv import load_dotenv

# Load environment variables from .env file, potentially in parent directory
load_dotenv("../../../.env") # Ensure this path is correct for your .env file

GEMINI_MODEL = os.getenv("GEMINI_MODEL", "gemini-2.5-flash-preview-05-20") # Load from env
instructions = read_markdown_file(ADMET_INSTRUCTIONS_MD_PATH)
instructions = f"""{instructions}"""


# Environment variables for ADMET MCP Server
# These should match the values used when running your ADMET_mcp_server.py
MCP_HOST = os.getenv("ADMET_HOST", "0.0.0.0") # Default from ADMET_mcp_server.py
MCP_PORT = os.getenv("ADMET_PORT", "8055")    # Default from ADMET_mcp_server.py
MCP_URL = f"http://{MCP_HOST}:{MCP_PORT}/sse"  # URL for the ADMET MCP server

root_agent = Agent(
    name="ADMETAgent",
    model=GEMINI_MODEL,
    instruction=instructions,
    tools=[
        MCPToolset(connection_params=SseServerParams(url=MCP_URL))
    ],
    description="Agent for Absorption, Distribution, Metabolism, Excretion, and Toxicity (ADMET) predictions using MCP server",
    after_tool_callback=after_tool_output_limit_callback)