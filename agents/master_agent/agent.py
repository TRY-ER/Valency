from custom_utils import MASTER_AGENT_INSTRUCTIONS_MD_PATH 
# Assuming 'custom_utils' is in the PYTHONPATH
from custom_utils.file_reader import read_markdown_file
from sub_agents import (
    protein_agent,
    drug_agent,
    drug_optimization_agent,
)
import os  # Added import
from google.adk.agents import Agent
from dotenv import load_dotenv
load_dotenv("../.env")

GEMINI_MODEL = os.getenv("GEMINI_MODEL", "gemini-2.5-flash-preview-05-20")  # Load from env
instructions = read_markdown_file(MASTER_AGENT_INSTRUCTIONS_MD_PATH)
instructions = f"""{instructions}"""

root_agent = Agent(
    name="MasterAgent",
    model=GEMINI_MODEL,
    instruction=instructions,
    sub_agents=[
        protein_agent,
        drug_agent,
        drug_optimization_agent,
    ],
    description="This agent is responsible for coordinating tasks across multiple sub-agents, including protein-related tasks, drug-related tasks, and drug optimization tasks."
)