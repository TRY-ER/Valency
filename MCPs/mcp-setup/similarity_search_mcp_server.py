import os
import base64
import io
from typing import List, Dict, Any, Optional
import json 
import sys 

from pydantic import BaseModel, Field 
import requests
import chromadb
from sentence_transformers import SentenceTransformer 
from rdkit import Chem
from rdkit.Chem import Draw 
from io import BytesIO 
from env_loader import load_env_vars 

from fastmcp import FastMCP # Changed import

load_env_vars() 

# --- Configuration ---
CHROMADB_PERSISTENT_PATH = os.getenv("CHROMADB_PERSISTENT_PATH", "./resources/chroma_store")
CHROMADB_SMILES_DB_NAME = os.getenv("CHROMADB_SMILES_DB_NAME", "smiles_data")
CHROMADB_PSMILES_DB_NAME = os.getenv("CHROMADB_PSMILES_DB_NAME", "psmiles_data")

SIMILARITY_MCP_HOST = os.getenv("SIMILARITY_MCP_HOST", "0.0.0.0")
SIMILARITY_MCP_PORT = int(os.getenv("SIMILARITY_MCP_PORT", "8059"))

ALLOWED_TYPES = ["MOL", "POLY", "PROT"]

# --- Pydantic Models (can be used for internal structuring) ---
class SimilarityQuery(BaseModel):
    input_type: str = Field(..., description="Type of input: MOL, POLY, or PROT")
    data: str = Field(..., description="Input data: SMILES string for MOL/POLY, PDB ID for PROT")
    k: int = Field(..., gt=0, description="Number of similar items to return")

class SearchResultItem(BaseModel):
    identifier: str
    score: float
    image: Optional[str] = None # Base64 encoded image string

class SimilarityResponse(BaseModel):
    status: str
    results: Optional[List[SearchResultItem]] = None
    message: Optional[str] = None

# --- Helper: ChromaDB Searcher ---
class ChromaSearcherReplicated:
    def __init__(self, collection_name: str, persist_directory: str, model_name: str = 'kuelumbus/polyBERT'):
        try:
            self.persist_directory = persist_directory
            self.client = chromadb.PersistentClient(path=self.persist_directory)
            # The collection will store vectors. Embeddings are generated before adding/querying.
            self.collection = self.client.get_or_create_collection(
                name=collection_name
                # No embedding_function needed here if we provide embeddings directly
            )
            self.model = SentenceTransformer(model_name)
            print(f"ChromaSearcherReplicated initialized for collection '{collection_name}' with model '{model_name}\'")
        except Exception as e:
            print(f"Error initializing ChromaSearcherReplicated for {collection_name}: {e}")
            raise RuntimeError(f"Failed to initialize ChromaDB or SentenceTransformer for {collection_name}: {str(e)}")

    def get_embedding(self, smiles: str) -> List[float]:
        try:
            return self.model.encode(smiles).tolist()
        except Exception as e:
            print(f"Error generating embedding for SMILES '{smiles}': {e}")
            raise RuntimeError(f"Failed to generate embedding: {str(e)}")

    def query(self, query_smiles: str, top_k: int) -> Dict[str, Any]:
        try:
            embedding = self.get_embedding(query_smiles)
            results = self.collection.query(
                query_embeddings=[embedding], # Query with the generated embedding
                n_results=top_k,
                include=['metadatas', 'distances'] # Ensure this matches expectations downstream
            )
            return results
        except Exception as e:
            print(f"Error querying ChromaDB collection with SMILES '{query_smiles}': {e}")
            raise RuntimeError(f"Error during similarity search in vector DB: {str(e)}")

# --- Helper: Molecule/Polymer Visualizer (RDKit) ---

# Define a base Visulizer class as it's inherited by MolVisualizer and PolymerVisualizer
class Visulizer:
    def visualize(self, source_str: str, size: tuple = (600, 600)):
        raise NotImplementedError

class MolVisualizer(Visulizer):

    def __init__(self):
        self.type = "MOL"

    def get_image(self, source_str: str, size: tuple = (600, 600)):
        mol = Chem.MolFromSmiles(source_str) # Correctly placed
        if not mol:
            print(f"Warning: Could not create RDKit molecule from SMILES: {source_str}")
            return None
        img = Draw.MolToImage(mol, size=size) 
        buffer = BytesIO()
        img.save(buffer, format="PNG")
        buffer.seek(0)
        # Ensure the output is a string for JSON serialization and correct indentation
        return f"data:image/png;base64,{base64.b64encode(buffer.read()).decode('utf-8')}"


    def visualize(self, source_str: str, size: tuple = (600, 600)):
        base64_image = self.get_image(source_str, size=size)
        return base64_image


class PolymerVisualizer(Visulizer):

    def __init__(self):
        self.type = "MOL" 

    def get_image(self, source_str: str, size: tuple = (600, 600)):
        mol = Chem.MolFromSmiles(source_str) # Correctly placed
        if not mol:
            print(f"Warning: Could not create RDKit polymer molecule from SMILES: {source_str}")
            return None
        img = Draw.MolToImage(mol, size=size)
        buffer = BytesIO()
        img.save(buffer, format="PNG")
        buffer.seek(0)
        # Ensure the output is a string for JSON serialization and correct indentation
        return f"data:image/png;base64,{base64.b64encode(buffer.read()).decode('utf-8')}"

    def visualize(self, source_str: str, size: tuple = (600, 600)):
        base64_image = self.get_image(source_str, size=size)
        return base64_image

# --- Helper: RCSB PDB Searcher & Image Fetcher ---
class RCSBQueryReplicated(BaseModel):
    entry_id: str
    rows: int

class RCSBSearcherReplicated:
    def search(self, query: RCSBQueryReplicated) -> Dict[str, Any]:
        # This is a simplified example. RCSB search can be complex.
        # This search finds other entries that were deposited with the same Uniprot Accession ID
        # as the primary Uniprot ID associated with the input PDB entry.
        # It's a form of relatedness, not necessarily structural similarity score.
        # For true structural similarity, a more complex query or different API endpoint is needed.
        
        search_url = "https://search.rcsb.org/rcsbsearch/v2/query"
        uniprot_query_payload = {
            "query": {
                "type": "terminal",
                "service": "text",
                "parameters": {"attribute": "rcsb_entry_container_identifiers.entry_id", "operator": "exact_match", "value": query.entry_id}
            },
            "request_options": {"return_all_hits": True},
            "return_type": "entry"
        }
        try:
            response_entry = requests.post(search_url, json=uniprot_query_payload, timeout=10)
            response_entry.raise_for_status()
            entry_data = response_entry.json()

            if not entry_data or entry_data.get("result_set") is None or not entry_data["result_set"]:
                 return {"result_set": []}

            uniprot_accession = None
            if entry_data["result_set"]:
                first_hit = entry_data["result_set"][0]
                if "rcsb_polymer_entity_container_identifiers" in first_hit and \
                   first_hit["rcsb_polymer_entity_container_identifiers"] and \
                   "reference_sequence_identifiers" in first_hit["rcsb_polymer_entity_container_identifiers"][0] and \
                   first_hit["rcsb_polymer_entity_container_identifiers"][0]["reference_sequence_identifiers"]:
                    for ref_id in first_hit["rcsb_polymer_entity_container_identifiers"][0]["reference_sequence_identifiers"]:
                        if ref_id.get("database_name") == "UniProt":
                            uniprot_accession = ref_id.get("database_accession")
                            break
            
            if not uniprot_accession:
                if query.rows >=1:
                     return {"result_set": [{"identifier": query.entry_id, "score": 1.0}]}
                return {"result_set": []}

            similarity_query_payload = {
                "query": {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                        "operator": "exact_match",
                        "value": uniprot_accession
                    }
                },
                "request_options": {
                    "pager": {"start": 0, "rows": query.rows},
                    "sort": [{"sort_by": "score", "direction": "desc"}]
                },
                "return_type": "entry"
            }
            
            response = requests.post(search_url, json=similarity_query_payload, timeout=15)
            response.raise_for_status()
            results = response.json()
            
            formatted_results = []
            if "result_set" in results:
                for item in results["result_set"]:
                    formatted_results.append({
                        "identifier": item["identifier"],
                        "score": item.get("score", 1.0) 
                    })
            return {"result_set": formatted_results}

        except requests.exceptions.RequestException as e:
            print(f"RCSB API request error: {e}")
            raise RuntimeError(f"Error communicating with RCSB PDB API: {str(e)}")
        except Exception as e:
            print(f"Error processing RCSB search: {e}")
            raise RuntimeError(f"Internal error during protein similarity search: {str(e)}")


def get_protein_image_replicated(pdb_id: str) -> Optional[str]:
    if not pdb_id or len(pdb_id) != 4: # Basic validation
        return None
    url = f"https://cdn.rcsb.org/images/structures/{pdb_id.lower()}_assembly-1.jpeg"
    try:
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            image_value = response.content
            if image_value:
                return f"data:image/jpeg;base64,{base64.b64encode(image_value).decode('utf-8')}"
        return None
    except requests.exceptions.RequestException as e:
        print(f"Failed to get image for PDB ID {pdb_id} from {url}: {e}")
        return None

# --- MCP Server Setup ---
mcp = FastMCP(
    "SimilaritySearchMCP",
    description="MCP Server for Molecule, Polymer, and Protein Similarity Searches using local logic."
    # Removed dependencies argument
    # host and port settings are not used here, passed to run() for SSE
)


def format_smiles_search_results_replicated(results: Dict[str, Any], input_type: str, k: int) -> List[Dict[str, Any]]: # Return List of Dicts
    formatted_items = []
    # ... existing logic ...
    if not results or not results.get("metadatas") or not results.get("distances"):
        return formatted_items

    metadatas_list = results["metadatas"][0] if results["metadatas"] and len(results["metadatas"]) > 0 else []
    distances_list = results["distances"][0] if results["distances"] and len(results["distances"]) > 0 else []
    
    num_results = min(len(metadatas_list), len(distances_list), k)

    for i in range(num_results):
        metadata = metadatas_list[i]
        distance = distances_list[i]
        
        smiles_key_str = "smiles" if "smiles" in metadata else "SMILES"
        identifier_smiles = metadata.get(smiles_key_str)
        
        if not identifier_smiles:
            continue

        # image_data = None
        # # Using the new visualizer classes
        # if input_type == "MOL":
        #     image_data = mol_viz.visualize(identifier_smiles, size=(300,300)) # Example size
        # elif input_type == "POLY":
        #     image_data = poly_viz.visualize(identifier_smiles, size=(300,300)) # Example size
        
        # Directly create dicts instead of SearchResultItem instances for MCP tool output
        formatted_items.append({
            "identifier": identifier_smiles,
            "score": 1.0 - distance,
            # "image": image_data
        })
    return formatted_items


def get_smiles_search_replicated_for_mcp(collection_name: str, query_smiles_data: str, k_val: int, input_type: str) -> Dict[str, Any]:
    """Helper function to perform SMILES search using ChromaSearcherReplicated."""
    try:
        searcher = ChromaSearcherReplicated(
            collection_name=collection_name, 
            persist_directory=CHROMADB_PERSISTENT_PATH
        )
        chroma_results = searcher.query(query_smiles=query_smiles_data, top_k=k_val)
        
        formatted_results = format_smiles_search_results_replicated(chroma_results, input_type, k_val)
        return {"status": "success", "data": formatted_results}
    except RuntimeError as e:
        # Errors from ChromaSearcherReplicated or embedding generation
        print(f"RuntimeError in SMILES search ({input_type}, {collection_name}): {e}")
        return {"status": "failed", "error": str(e)}
    except Exception as e:
        # Catch any other unexpected errors
        print(f"Unexpected error in get_smiles_search_replicated_for_mcp ({input_type}, {collection_name}): {e}")
        return {"status": "failed", "error": f"An unexpected error occurred during {input_type} similarity search: {str(e)}"}


@mcp.tool()
def perform_similarity_search(input_type: str, data: str, k: int) -> str:
    """
    Performs a similarity search for Molecules (MOL), Polymers (POLY), or Proteins (PROT).
    For MOL/POLY, uses SMILES strings and a local vector DB.
    For PROT, uses PDB IDs and queries RCSB PDB.

    Args:
        input_type: Type of input: "MOL", "POLY", or "PROT".
        data: Input data: SMILES string for MOL/POLY, PDB ID for PROT.
        k: Number of similar items to return (must be > 0).

    Returns:
        A JSON string containing the search results or an error message.
        Success: {"status": "success", "data": [{"identifier": "...", "score": 0.0, "image": "base64_string_or_null"}, ...]}
        Failure: {"status": "failed", "error": "Error message"}
    """
    if k < 1:
        return json.dumps({"status": "failed", "error": "k should be a positive integer"})
    if input_type not in ALLOWED_TYPES:
        return json.dumps({"status": "failed", "error": f"Input type should be one of {ALLOWED_TYPES}"})

    response_dict: Dict[str, Any] = {}

    try:
        if input_type == "MOL":
            if not CHROMADB_SMILES_DB_NAME:
                return json.dumps({"status": "failed", "error": "MOL collection name not configured."})
            # Pass 'data' (which is the SMILES string from the tool input) as query_smiles_data
            response_dict = get_smiles_search_replicated_for_mcp(CHROMADB_SMILES_DB_NAME, data, k, "MOL")
        
        elif input_type == "POLY":
            if not CHROMADB_PSMILES_DB_NAME:
                return json.dumps({"status": "failed", "error": "POLY collection name not configured."})
            # Pass 'data' as query_smiles_data
            response_dict = get_smiles_search_replicated_for_mcp(CHROMADB_PSMILES_DB_NAME, data, k, "POLY")
        
        elif input_type == "PROT":
            rcsb_query = RCSBQueryReplicated(entry_id=data, rows=k) # Use the Pydantic model internally
            searcher = RCSBSearcherReplicated()
            rcsb_results_dict = searcher.search(rcsb_query)
            
            returnable_items: List[Dict[str, Any]] = []
            if "result_set" in rcsb_results_dict:
                for r_item in rcsb_results_dict["result_set"][:k]:
                    pdb_id = r_item.get("identifier")
                    if not pdb_id:
                        continue
                    # image_data = get_protein_image_replicated(pdb_id)
                    returnable_items.append({
                        "identifier": pdb_id,
                        "score": r_item.get("score", 0.0),
                        # "image": image_data
                    })
            response_dict = {"status": "success", "data": returnable_items}

        else: # Should be caught by initial validation
            response_dict = {"status": "failed", "error": "Invalid input_type specified."}

    except RuntimeError as e: # Catch errors raised by helpers
        print(f"RuntimeError in perform_similarity_search tool: {e}")
        response_dict = {"status": "failed", "error": str(e)}
    except requests.exceptions.RequestException as e: # Specifically for RCSB network issues if not caught by helper
        print(f"RequestException in PROT similarity search: {e}")
        response_dict = {"status": "failed", "error": f"RCSB API Communication Error: {str(e)}"}
    except Exception as e:
        print(f"Unexpected error in perform_similarity_search tool: {e}")
        response_dict = {"status": "failed", "error": f"An unexpected error occurred: {str(e)}"}
        
    return json.dumps(response_dict)

# --- Main execution ---
if __name__ == "__main__":
    print(f"Starting Similarity Search MCP Server...")
    print(f"Allowed input types: {ALLOWED_TYPES}")
    print(f"ChromaDB Path: {CHROMADB_PERSISTENT_PATH}")
    print(f"ChromaDB MOL Collection: {CHROMADB_SMILES_DB_NAME}")
    print(f"ChromaDB POLY Collection: {CHROMADB_PSMILES_DB_NAME}")
    
    if CHROMADB_PERSISTENT_PATH != ":memory:":
        os.makedirs(CHROMADB_PERSISTENT_PATH, exist_ok=True)
    
    # transport = os.getenv("MCP_TRANSPORT", "sse") # No longer needed for direct SSE run
    print(f"MCP Server Name: {mcp.name}")
    sys.stdout.flush()
    
    # Default to SSE transport with host and port passed directly
    print(f"Attempting to run Similarity Search MCP Server with FastMCP SSE transport on host {SIMILARITY_MCP_HOST}, port {SIMILARITY_MCP_PORT}")
    sys.stdout.flush()
    mcp.run(transport="sse", host=SIMILARITY_MCP_HOST, port=SIMILARITY_MCP_PORT)
    # else:
    #     # Fallback or error for unknown transport
    #     print(f"Unknown transport: {transport}. Supported: 'stdio', 'sse'. Defaulting to stdio.")
    #     sys.stdout.flush()
    #     mcp.run(transport="stdio")

