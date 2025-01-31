from fastapi import APIRouter
from engine.explorers.mol_explorer.mol_info import MoleculeInfo 
from engine.explorers.prot_explorer.prot_info import ProteinInfo 
from engine.explorers.poly_explorer.poly_info import PolymerInfo
from engine.base import ALLOWED_TYPES

router = APIRouter()

informer_mapper = {
    "MOL": MoleculeInfo(),
    "PROT": ProteinInfo(),
    "POLY": PolymerInfo()
}

@router.post("/")
def get_molecule_info(data: dict):
    informer_type = data["type"]
    if informer_type is None or informer_type == "" or informer_type not in ALLOWED_TYPES:
        return {"status": "error",
                "message": f"Visualization type should be one of {ALLOWED_TYPES}"}
    
    informer = informer_mapper.get(informer_type)
    return {
        "status": "success",
        "info": informer.get_info(str(data["value"]))
    }