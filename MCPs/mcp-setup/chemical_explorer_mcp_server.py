import os
import json
import sys 
from env_loader import load_env_vars

load_env_vars()

from fastmcp import FastMCP # Changed import

# Attempt to import RDKit
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    print("Warning: RDKit not found. Chemical property calculations will not be available.")
    # Define dummy Descriptors and Chem if RDKit is not available to avoid runtime errors on class definition
    class Descriptors:
        @staticmethod
        def MolWt(mol): return 0.0  # Return float
        @staticmethod
        def RingCount(mol): return 0  # Return int
        @staticmethod
        def HeavyAtomCount(mol): return 0  # Return int
        @staticmethod
        def NumHDonors(mol): return 0  # Return int
        @staticmethod
        def NumHAcceptors(mol): return 0  # Return int
        @staticmethod
        def NumRotatableBonds(mol): return 0  # Return int
        @staticmethod
        def TPSA(mol): return 0.0  # Return float

    class Chem:
        @staticmethod
        def MolFromSmiles(smiles_string): return None # Remains None, as this is checked
        @staticmethod
        def rdMolDescriptors_CalcMolFormula(mol): return "" # Return empty string
        class rdMolDescriptors: # Nested class needs to be defined
            @staticmethod
            def CalcMolFormula(mol): return "" # Return empty string

# Attempt to import pypdb
try:
    from pypdb import get_info as get_pdb_info
    PYPDB_AVAILABLE = True
except ImportError:
    PYPDB_AVAILABLE = False
    print("Warning: pypdb not found. Protein information retrieval will not be available.")
    # Define a dummy get_pdb_info if pypdb is not available
    def get_pdb_info(pdb_id_str: str): return None


# --- Environment Variables for Chemical Explorer MCP Server ---
chemical_explorer_mcp_host = os.getenv("CHEMICAL_EXPLORER_MCP_HOST", "0.0.0.0")
chemical_explorer_mcp_port = int(os.getenv("CHEMICAL_EXPLORER_MCP_PORT", "8055"))

mcp = FastMCP(
    "Chemical Explorer MCP Server",
    # Removed dependencies from here, ensure they are installed in the environment if needed by tools
    # Removed host/port from settings, will be passed to run()
    description="MCP server for retrieving information about molecules, polymers using RDKit and proteins using PDB."
)

# --- Self-contained Classes for Molecule and Polymer using RDKit ---

class MoleculeInfoProvider:
    def __init__(self, smiles: str):
        self.smiles = smiles
        self.mol = Chem.MolFromSmiles(smiles) if RDKIT_AVAILABLE else None

    def get_info(self) -> dict:
        if not RDKIT_AVAILABLE:
            return {"error": "RDKit is not installed. Cannot provide molecule information."}
        if not self.mol:
            return {"error": f"Invalid SMILES string: '{self.smiles}'"}

        info = {
            "smiles": self.smiles,
            "Molecular Formula": Chem.rdMolDescriptors.CalcMolFormula(self.mol),
            "Molecular Weight": round(float(Descriptors.MolWt(self.mol)), 3),
            "Heavy Atoms Count": Descriptors.HeavyAtomCount(self.mol),
            "H Bond Donor Count": Descriptors.NumHDonors(self.mol),
            "H Bond Acceptor Count": Descriptors.NumHAcceptors(self.mol),
            "Rotatable Bonds Count": Descriptors.NumRotatableBonds(self.mol),
            "TPSA": round(Descriptors.TPSA(self.mol), 2),
            "Number of Rings": Descriptors.RingCount(self.mol),
        }
        return info

class PolymerInfoProvider:
    def __init__(self, psmiles: str):
        self.psmiles = psmiles
        self.mol = Chem.MolFromSmiles(psmiles) if RDKIT_AVAILABLE else None

    def get_info(self) -> dict:
        if not RDKIT_AVAILABLE:
            return {"error": "RDKit is not installed. Cannot provide polymer information."}
        if not self.mol:
            return {"error": f"Invalid Polymer SMILES (PSMILES) string: '{self.psmiles}'"}

        # Get indexes of [*] content in the PSMILES string
        wildcard_indices = []
        if self.mol: # Ensure mol is not None before iterating
            wildcard_indices = [atom.GetIdx() for atom in self.mol.GetAtoms() if atom.GetSymbol() == '*']
        wildcard_indices_str = ",".join([str(index) for index in wildcard_indices])

        info = {
            "psmiles": self.psmiles,
            "Monomer Molecular Formula": Chem.rdMolDescriptors.CalcMolFormula(self.mol),
            "Monomer Molecular Weight": round(float(Descriptors.MolWt(self.mol)), 3),
            "Number of Rings in Monomer": Descriptors.RingCount(self.mol),
            "Open Bond Indexes": wildcard_indices_str
        }
        return info

# --- Self-contained Class for Protein Information using pypdb ---
class PDBProteinInfoProvider:
    def __init__(self, pdb_id: str):
        self.pdb_id = pdb_id

    def _extract_simplified_pdb_data(self, pdb_data: dict) -> dict:
        """
        Extracts significant scientific properties from a complex PDB data dictionary.
        """
        if not pdb_data: # Handle case where pdb_data might be None (e.g. pypdb not found or invalid ID)
            return {"error": f"Could not retrieve data for PDB ID: '{self.pdb_id}'"}

        simplified_data = {}
        simplified_data['Pdb_Id'] = pdb_data.get('rcsb_id') or pdb_data.get('entry', {}).get('id')
        simplified_data['Title'] = pdb_data.get('struct', {}).get('title')
        authors_list = [author['name'] for author in pdb_data.get('audit_author', [])]
        simplified_data['Authors'] = ', '.join(authors_list) if authors_list else None
        citation = pdb_data.get('citation', [])
        if citation:
            primary_citation = next((c for c in citation if c.get('id') == 'primary'), citation[0])
            simplified_data['Journal'] = primary_citation.get('rcsb_journal_abbrev') or primary_citation.get('journal_abbrev')
            simplified_data['Year'] = primary_citation.get('year')
            simplified_data['Volume'] = primary_citation.get('journal_volume')
            page_first = primary_citation.get('page_first')
            page_last = primary_citation.get('page_last')
            if page_first and page_last:
                pages = f"{page_first}-{page_last}"
            elif page_first:
                pages = page_first
            else:
                pages = None
            simplified_data['Pages'] = pages
            simplified_data['Doi'] = primary_citation.get('pdbx_database_id_doi')
            simplified_data['Pubmed_Id'] = primary_citation.get('pdbx_database_id_pub_med')
        else:
            simplified_data['Journal'] = None
            simplified_data['Year'] = None
            simplified_data['Volume'] = None
            simplified_data['Pages'] = None
            simplified_data['Doi'] = None
            simplified_data['Pubmed_Id'] = None
        exptl = pdb_data.get('exptl', [])
        if exptl:
            simplified_data['Experiment_Method'] = exptl[0].get('method')
        else:
            simplified_data['Experiment_Method'] = None
        simplified_data['Molecular_Weight_(kDa)'] = pdb_data.get('rcsb_entry_info', {}).get('molecular_weight')
        simplified_data['Deposited_Model_Count'] = pdb_data.get('rcsb_entry_info', {}).get('deposited_model_count')
        simplified_data['Polymer_entity_count'] = pdb_data.get('rcsb_entry_info', {}).get('polymer_entity_count')
        simplified_data['Polymer_monomer_count'] = pdb_data.get('rcsb_entry_info', {}).get('deposited_polymer_monomer_count')
        simplified_data['Structural_Features'] = pdb_data.get('struct_keywords', {}).get('text')
        release_date = pdb_data.get('rcsb_accession_info', {}).get('initial_release_date', '')
        if 'T' in release_date:
            simplified_data['Release_Date'] = release_date.split('T')[0]
        else:
            simplified_data['Release_Date'] = release_date
        resolution = pdb_data.get('rcsb_entry_info', {}).get('resolution_combined', [None])
        if resolution and resolution[0] is not None:
            simplified_data['Resolution'] = resolution[0]
        else:
            simplified_data['Resolution'] = None
        return simplified_data

    def get_info(self) -> dict:
        if not PYPDB_AVAILABLE:
            return {"error": "pypdb is not installed. Cannot provide protein information."}
        try:
            raw_pdb_data = get_pdb_info(self.pdb_id)
            if not raw_pdb_data: # Check if pypdb returned None (e.g. invalid PDB ID)
                 return {"error": f"No data found for PDB ID: '{self.pdb_id}'. It might be an invalid ID."}
            return self._extract_simplified_pdb_data(raw_pdb_data)
        except Exception as e:
            return {"error": f"Failed to retrieve or process PDB data for '{self.pdb_id}': {str(e)}"}


# --- Chemical Explorer API Tools ---

@mcp.tool()
def get_molecule_information(smiles_string: str) -> str:
    """
    Retrieves RDKit-calculated information about a specific molecule from its SMILES string.
    Args:
        smiles_string: The SMILES string of the molecule (e.g., 'CCO', 'c1ccccc1C(=O)O').
    Returns:
        A JSON string containing the molecule's details, or an error message.
    """
    if not RDKIT_AVAILABLE:
        return json.dumps({"error": "RDKit backend not available on the server."})
    
    provider = MoleculeInfoProvider(smiles_string)
    info = provider.get_info()
    return json.dumps(info)

@mcp.tool()
def get_polymer_information(psmiles_string: str) -> str:
    """
    Retrieves RDKit-calculated information about a specific polymer repeating unit from its PSMILES string.
    Args:
        psmiles_string: The Polymer SMILES (PSMILES) string of the repeating unit (e.g., '[*]CC[*]', '[*]NCCCCN[*]C(=O)CCCC(=O)').
    Returns:
        A JSON string containing the polymer unit's details, or an error message.
    """
    if not RDKIT_AVAILABLE:
        return json.dumps({"error": "RDKit backend not available on the server."})

    provider = PolymerInfoProvider(psmiles_string)
    info = provider.get_info()
    return json.dumps(info)

@mcp.tool()
def get_protein_information(pdb_id_string: str) -> str:
    """
    Retrieves information about a specific protein from its PDB ID using pypdb.
    Args:
        pdb_id_string: The PDB ID of the protein (e.g., '6M0J', '1TIM').
    Returns:
        A JSON string containing the protein's details, or an error message.
    """
    if not PYPDB_AVAILABLE:
        return json.dumps({"error": "pypdb backend not available on the server."})
    
    provider = PDBProteinInfoProvider(pdb_id_string)
    info = provider.get_info()
    return json.dumps(info)

# --- Main execution for direct run or mcp dev ---
if __name__ == "__main__":
    if not RDKIT_AVAILABLE:
        print("CRITICAL: RDKit library is not installed. This MCP server requires RDKit to function for molecule/polymer tools.")
        print("Please install it, e.g., using: pip install rdkit")
        sys.stdout.flush()
        # Optionally, exit if RDKit is critical and not found for certain core functionalities
        # exit(1)
    
    if not PYPDB_AVAILABLE:
        print("WARNING: pypdb library is not installed. Protein information retrieval will not be available.")
        print("Please install it, e.g., using: pip install pypdb")
        sys.stdout.flush()

    # transport = os.getenv("MCP_TRANSPORT", "sse") # No longer needed for direct SSE run
    print(f"Starting Chemical Explorer MCP Server...")
    print(f"Server Name: {mcp.name}")
    sys.stdout.flush() 

    # if transport == "stdio":
    #     print("Running Chemical Explorer MCP Server with stdio transport")
    #     sys.stdout.flush() 
    #     mcp.run(transport="stdio") # Still allow stdio if explicitly set via env var for testing
    # elif transport == "sse":
    # Default to SSE transport with host and port passed directly
    print(f"Attempting to run Chemical Explorer MCP Server with FastMCP SSE transport on host {chemical_explorer_mcp_host}, port {chemical_explorer_mcp_port}")
    sys.stdout.flush() 
    mcp.run(transport="sse", host=chemical_explorer_mcp_host, port=chemical_explorer_mcp_port)
    # else:
    #     raise ValueError(f"Unknown transport: {transport}. Supported: 'stdio', 'sse'.")

