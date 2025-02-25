from engine.explorers.base import Visulizer
from rdkit import Chem
from rdkit.Chem import Draw
from io import BytesIO
import base64


class MolVisualizer(Visulizer):

    def __init__(self):
        self.type = "MOL"

    def get_image(self, source_str: str, size: tuple = (600, 600)):
        mol = Chem.MolFromSmiles(source_str)
        # mol = Chem.AddHs(mol)
        if not mol:
            return None
        img = Draw.MolToImage(mol, size=(600, 600))
        buffer = BytesIO()
        img.save(buffer, format="PNG")
        buffer.seek(0)
        return base64.b64encode(buffer.read()).decode("utf-8")

    def visualize(self, source_str: str, size: tuple = (600, 600)):
        base64_image = self.get_image(source_str, size=size)
        return base64_image
