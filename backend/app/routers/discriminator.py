from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from engine.vector_stores.smiles_vectorstore_chromadb import ChromaSearcher
from get_env_vars import CHROMADB_PERSISTENT_PATH, CHROMADB_SMILES_DB_NAME, CHROMADB_PSMILES_DB_NAME
from engine.base import ALLOWED_TYPES
from engine.discriminator.similarity_w_rcsb import RCSBSearcher, RCSBQuery
from requests.exceptions import HTTPError
from engine.explorers.mol_explorer.mol_visualizer import MolVisualizer
from engine.explorers.poly_explorer.poly_visualizer import PolymerVisualizer 
import requests
import base64

router = APIRouter()

class SimilarityQuery(BaseModel):
    input_type: str
    data: str
    k: int

def format_smiles_search_results(results, input_type: str):
    formatted_results = []
    # try:
    smiles_key_str = "smiles" if "smiles" in results["metadatas"][0][0] else "SMILES" 
    for metadata, distance in zip(results["metadatas"][0], results["distances"][0]):
        if input_type == "MOL":
            image_data = MolVisualizer().visualize(metadata.get(smiles_key_str), (100, 100))
        elif input_type == "POLY":
            image_data = PolymerVisualizer().visualize(metadata.get(smiles_key_str), (100, 100))
        else:
            raise HTTPException(status_code=400, detail=f"Validation type should be one of {ALLOWED_TYPES}")
        formatted_results.append({
            "identifier": metadata.get(smiles_key_str),
            "score": distance,
            "image": image_data
        })
    return {"status": "success", "results": formatted_results} 

def get_smiles_search(collection_name, payload, input_type: str):
    searcher = ChromaSearcher(collection_name=collection_name, persist_directory=CHROMADB_PERSISTENT_PATH)
    results = searcher.query(payload.data, top_k=payload.k)
    return format_smiles_search_results(results, input_type)

def get_protein_image(pdb_id: str):
    url = f"https://cdn.rcsb.org/images/structures/{pdb_id.lower()}_assembly-1.jpeg"
    # get the image using url and return the image data as base64 encoding

    image_value = None 
    response = requests.get(url)
    if response.status_code == 200:
        image_value = response.content
        if image_value is None:
            raise HTTPException(status_code=400, detail="Failed to get image data")
        image_base64 = base64.b64encode(image_value).decode("utf-8")
        return image_base64


@router.post("/ssearch/local")
def similarity_search(payload: SimilarityQuery):
    if payload.k < 1:
        raise HTTPException(status_code=400, detail="k should be a positive integer")
    if payload.input_type not in ALLOWED_TYPES:
        raise HTTPException(status_code=400, detail=f"Validation type should be one of {ALLOWED_TYPES}")
    if payload.input_type == "MOL":
        collection_name = CHROMADB_SMILES_DB_NAME 
        return get_smiles_search(collection_name, payload, "MOL")
    elif payload.input_type == "POLY":
        collection_name = CHROMADB_PSMILES_DB_NAME
        return get_smiles_search(collection_name, payload, "POLY")
    elif payload.input_type == "PROT":
        query = RCSBQuery(entry_id=payload.data, rows=payload.k)
        searcher = RCSBSearcher()
        try:
            results = searcher.search(query)
            returnable = [] 
            for r in results["result_set"]:
                image_data = get_protein_image(r["identifier"])
                returnable.append({
                    "identifier": r["identifier"],
                    "score": r["score"],
                    "image": image_data
                })
            return {"status": "success", "results": returnable}
        except HTTPError as e: 
            return {"status": "failed", "message": "Try with a valid PDB Id"} 
    