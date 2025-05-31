# filepath: /home/kalki/src/valency/agents/drug_optimization_agent/agent.py
from custom_utils.repo.declare import DRUG_OPTIMIZATION_INSTRUCTIONS_MD_PATH
from custom_utils.file_reader import read_markdown_file
from ..sub_sub_agents import (
    admet_agent,
    brics_agent,
)
import os  # Added import
from google.adk.agents import Agent
from dotenv import load_dotenv

load_dotenv("../../.env")  # Adjusted path for .env

GEMINI_MODEL = os.getenv("GEMINI_MODEL", "gemini-2.5-flash-preview-05-20")  # Load from env
instructions = read_markdown_file(DRUG_OPTIMIZATION_INSTRUCTIONS_MD_PATH)
instructions = f"""{instructions}"""

root_agent = Agent(
    name="DrugOptimizationAgent",
    model=GEMINI_MODEL,
    instruction=instructions,
    sub_agents=[
        admet_agent,
        brics_agent
    ],
    description="This agent is responsible for drug optimization tasks, including ADMET predictions and BRICS decomposition."
)
