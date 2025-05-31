from custom_utils import PUBCHEM_INSTRUCTIONS_MD_PATH 
import asyncio
import os
from google.adk.tools.mcp_tool.mcp_toolset import (
    MCPToolset,
    SseServerParams
)

# Assuming 'custom_utils' is in the PYTHONPATH
from custom_utils.file_reader import read_markdown_file
from google.adk.agents import Agent
from dotenv import load_dotenv
load_dotenv("../.env")

GEMINI_MODEL = os.getenv("GEMINI_MODEL", "gemini-2.5-flash-preview-05-20")
instructions = read_markdown_file(PUBCHEM_INSTRUCTIONS_MD_PATH)
instructions = f"""{instructions}"""


# Default to localhost if not set
MCP_HOST = os.getenv("PUBCHEM_HOST", "localhost")
MCP_PORT = os.getenv("PUBCHEM_PORT", "8054")  # Default to 8058 if not set
MCP_URL = f"http://{MCP_HOST}:{MCP_PORT}/sse"  # URL for the MCP server

root_agent = Agent(
    name="AlphafoldAgent",
    model=GEMINI_MODEL,
    instruction=instructions,
    tools=[
        MCPToolset(connection_params=SseServerParams(url=MCP_URL))
    ])