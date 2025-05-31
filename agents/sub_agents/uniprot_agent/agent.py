from custom_utils import  UNIPROT_INSTRUCTIONS_MD_PATH# This will be a new constant
import os
from google.adk.tools.mcp_tool.mcp_toolset import (
    MCPToolset,
    SseServerParams
)
from google.adk.agents import Agent
from custom_utils.file_reader import read_markdown_file
from dotenv import load_dotenv

# Load environment variables from a .env file located one level above the current script's directory
# This is useful if your .env file is in the parent directory of the 'sub_agents' folder
load_dotenv(dotenv_path=os.path.join(os.path.dirname(__file__), "../../.env"))

GEMINI_MODEL = os.getenv("GEMINI_MODEL", "gemini-1.5-flash-latest") # Use a default if not set

instructions = read_markdown_file(UNIPROT_INSTRUCTIONS_MD_PATH)
instructions = f"""{instructions}"""


# Default to localhost and standard UniProt MCP port if not set in .env
# These environment variables should be set where your MCP server is running
MCP_HOST = os.getenv("UNIPROT_HOST", "localhost")
MCP_PORT = os.getenv("UNIPROT_PORT", "8052")  # Default UniProt MCP server port
MCP_URL = f"http://{MCP_HOST}:{MCP_PORT}/sse"  # URL for the MCP server

root_agent = Agent(
    name="UniprotAgent",
    model=GEMINI_MODEL,
    instruction=instructions,
    tools=[
        MCPToolset(connection_params=SseServerParams(url=MCP_URL))
    ],
    description="An agent for UniProt protein and gene understanding from UniProt ID or string query",)

# To make this agent runnable, you might want to add something like:
# if __name__ == '__main__':
#     # Example of how to run the agent (requires appropriate input handling)
#     # from google.adk.agents import run_agent_loop
#     # run_agent_loop(root_agent)
#     print(f"UniprotAgent initialized with MCP URL: {MCP_URL}")
#     print("To run this agent, you would typically integrate it into a larger application or use a runner.")
