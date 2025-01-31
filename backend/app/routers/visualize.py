from fastapi import APIRouter
from engine.explorers.mol_explorer.mol_visualizer import MolVisualizer
from engine.explorers.prot_explorer.prot_visualizer import ProtVisualizer
from engine.explorers.poly_explorer.poly_visualizer import PolymerVisualizer 
from engine.base import ALLOWED_TYPES

router = APIRouter()

visualizer_mapper = {
    "MOL": MolVisualizer(),
    "PROT": ProtVisualizer(),
    "POLY": PolymerVisualizer()
}

@router.post("/")
def visualize_molecule(data: dict):
    visualizer_type = data["type"]
    if visualizer_type is None or visualizer_type == "" or visualizer_type not in ALLOWED_TYPES:
        return {"status": "error",
                "message": f"Visualization type should be one of {ALLOWED_TYPES}"}
    visualizer = visualizer_mapper.get(visualizer_type)
    return {"status": "success",
            "image_data": visualizer.visualize(str(data["value"]))}