from custom_utils import BRICS_INSTRUCTIONS_MD_PATH
from custom_utils import after_tool_output_limit_callback
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
load_dotenv("../../../.env")

GEMINI_MODEL = os.getenv("GEMINI_MODEL", "gemini-2.5-flash-preview-05-20")
instructions = read_markdown_file(BRICS_INSTRUCTIONS_MD_PATH)
instructions = f"""{instructions}"""


# Default to localhost if not set
MCP_HOST = os.getenv("BRICS_HOST", "localhost")
MCP_PORT = os.getenv("BRICS_PORT", "8058")  # Default to 8058 if not set
MCP_URL = f"http://{MCP_HOST}:{MCP_PORT}/sse"  # URL for the MCP server

root_agent = Agent(
    name="BRICSAgent",
    model=GEMINI_MODEL,
    instruction=instructions,
    tools=[
        MCPToolset(connection_params=SseServerParams(url=MCP_URL))
    ],
    description="Agent for breaking of retrosynthetically interesting chemical substructures (BRICS) based generation using MCP server",
    after_tool_callback=after_tool_output_limit_callback)