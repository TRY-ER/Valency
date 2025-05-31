import os

# Get the absolute path of the directory where this script (declare.py) is located.
# This is /home/kalki/src/valency/agents/custom_utils/repo/
REPO_DIR_PATH = os.path.dirname(os.path.abspath(__file__))

# Declare variables for the absolute paths of markdown files in this directory.
# Example for sample_instructions.md:
SAMPLE_INSTRUCTIONS_MD_PATH = os.path.join(REPO_DIR_PATH, "sample_instructions.md")
RCSB_INSTRUCTIONS_MD_PATH = os.path.join(REPO_DIR_PATH, "rcsb_instructions.md")
ALPHAFOLD_INSTRUCTIONS_MD_PATH = os.path.join(REPO_DIR_PATH, "alphafold_instructions.md")
BRICS_INSTRUCTIONS_MD_PATH = os.path.join(REPO_DIR_PATH, "brics_instructions.md")
ADMET_INSTRUCTIONS_MD_PATH = os.path.join(REPO_DIR_PATH, "admet_instructions.md")
CHEMBL_INSTRUCTIONS_MD_PATH = os.path.join(REPO_DIR_PATH, "chembl_instructions.md")
PUBCHEM_INSTRUCTIONS_MD_PATH = os.path.join(REPO_DIR_PATH, "pubchem_instructions.md")
SS_INSTRUCTIONS_MD_PATH = os.path.join(REPO_DIR_PATH, "similar_search_instructions.md")
UNIPROT_INSTRUCTIONS_MD_PATH = os.path.join(REPO_DIR_PATH, "uniprot_instructions.md")

# nest agent instructions
PROTEIN_AGENT_INSTRUCTIONS_MD_PATH = os.path.join(REPO_DIR_PATH, "protein_agent_instructions.md")
DRUG_AGENT_INSTRUCTIONS_MD_PATH = os.path.join(REPO_DIR_PATH, "drug_agent_instructions.md")
DRUG_OPTIMIZATION_INSTRUCTIONS_MD_PATH  = os.path.join(REPO_DIR_PATH, "drug_optimization_instructions.md")

# master agent instructions
MASTER_AGENT_INSTRUCTIONS_MD_PATH = os.path.join(REPO_DIR_PATH, "master_agent_instructions.md")