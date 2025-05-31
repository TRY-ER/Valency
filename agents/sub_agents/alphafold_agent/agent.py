from custom_utils import ALPHAFOLD_INSTRUCTIONS_MD_PATH
import asyncio
import os
from google.adk.tools.mcp_tool.mcp_toolset import (
    MCPToolset,
    SseServerParams
)
from google.adk.tools import ToolContext
from google.adk.agents import LoopAgent

# Assuming 'custom_utils' is in the PYTHONPATH
from custom_utils.file_reader import read_markdown_file
from google.adk.agents import Agent
from dotenv import load_dotenv
load_dotenv("../.env")

GEMINI_MODEL = "gemini-2.0-flash"
instructions = read_markdown_file(ALPHAFOLD_INSTRUCTIONS_MD_PATH)
instructions = f"""{instructions}"""


# Default to localhost if not set
MCP_HOST = os.getenv("ALPHAFOLD_HOST", "localhost")
MCP_PORT = os.getenv("ALPHAFOLD_PORT", "8056")  # Default to 8058 if not set
MCP_URL = f"http://{MCP_HOST}:{MCP_PORT}/sse"  # URL for the MCP server

root_agent = Agent(
    name="AlphafoldAgent",
    model=GEMINI_MODEL,
    instruction=instructions,
    tools=[
        MCPToolset(connection_params=SseServerParams(url=MCP_URL))
    ],
    description="An agent for AlphaFold protein structure prediction from UniPROT key",)