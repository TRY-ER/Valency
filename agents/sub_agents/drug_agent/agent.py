from custom_utils import DRUG_AGENT_INSTRUCTIONS_MD_PATH 
# Assuming 'custom_utils' is in the PYTHONPATH
from custom_utils.file_reader import read_markdown_file
from ..sub_sub_agents import (
    pubchem_agent,
    chembl_agent,
)
import os
from google.adk.agents import Agent
from dotenv import load_dotenv
load_dotenv("../../.env")

GEMINI_MODEL = os.getenv("GEMINI_MODEL", "gemini-2.5-flash-preview-05-20")
instructions = read_markdown_file(DRUG_AGENT_INSTRUCTIONS_MD_PATH)
instructions = f"""{instructions}"""

root_agent = Agent(
    name="DrugAgent",
    model=GEMINI_MODEL,
    instruction=instructions,
    sub_agents=[
        pubchem_agent,
        chembl_agent
    ],
    description="This agent is responsible for drug-related tasks, including searching for drugs in PubChem and ChEMBL databases."
)