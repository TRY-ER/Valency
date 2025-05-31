from custom_utils import PROTEIN_AGENT_INSTRUCTIONS_MD_PATH
# Assuming 'custom_utils' is in the PYTHONPATH
from custom_utils.file_reader import read_markdown_file
from ..sub_sub_agents import (
    uniprot_agent,
    rcsb_pdb_agent,
    alphafold_agent
)
import os  # Added import
from google.adk.agents import Agent
from dotenv import load_dotenv
load_dotenv("../../.env")

GEMINI_MODEL = os.getenv("GEMINI_MODEL", "gemini-2.5-flash-preview-05-20")  # Load from env
instructions = read_markdown_file(PROTEIN_AGENT_INSTRUCTIONS_MD_PATH)
instructions = f"""{instructions}"""

root_agent = Agent(
    name="ProteinAgent",
    model=GEMINI_MODEL,
    instruction=instructions,
    sub_agents=[
        uniprot_agent,
        rcsb_pdb_agent,
        alphafold_agent
    ],
    description="This agent is responsible for protein-related tasks, including searching for proteins in UniProt, RCSB PDB, and AlphaFold databases."
)