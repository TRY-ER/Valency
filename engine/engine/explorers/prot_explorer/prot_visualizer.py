from engine.explorers.base import Visulizer
from Bio import SeqIO
from Bio.PDB import *
from Bio.SeqUtils import ProtParam
import matplotlib.pyplot as plt
from io import BytesIO
import base64

class ProtVisualizer(Visulizer):
    def __init__(self):
        self.type = "PROT"
    
    def get_image(self, source_str: str):
        try:
            # Create a ProteinAnalysis object
            protein = ProtParam.ProteinAnalysis(source_str)
            
            # Get secondary structure prediction
            secondary_struct = protein.secondary_structure_fraction()
            
            # Create a simple visualization using matplotlib
            labels = ['Helix', 'Turn', 'Sheet']
            sizes = [secondary_struct[0], secondary_struct[1], secondary_struct[2]]
            
            # Create pie chart
            plt.figure(figsize=(8, 8))
            plt.pie(sizes, labels=labels, autopct='%1.1f%%')
            plt.title('Predicted Secondary Structure Distribution')
            
            # Save to buffer
            buffer = BytesIO()
            plt.savefig(buffer, format='png')
            buffer.seek(0)
            plt.close()
            
            return base64.b64encode(buffer.getvalue()).decode('utf-8')
            
        except Exception as e:
            return None
    
    def visualize(self, source_str: str):
        base64_image = self.get_image(source_str)
        return base64_image