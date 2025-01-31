from rdkit import Chem
from rdkit.Chem import Descriptors

class MoleculeInfo:
    def get_info(self, source_str: str) -> dict:
        """
        Return basic RDKit descriptors from a SMILES string.
        """
        mol = Chem.MolFromSmiles(source_str)
        if not mol:
            return {"error": "Invalid SMILES"}

        # Basic properties about the molecule:
        # MolecularFormula  - chemical formula for the molecule
        # MolecularWeight   - approximate weight in daltons
        # NumHeavyAtoms     - count of non-hydrogen atoms
        # NumHBD            - hydrogen bond donor count
        # NumHBA            - hydrogen bond acceptor count
        # NumRotatableBonds - number of freely rotating bonds
        # TPSA              - topological polar surface area
        # NumRings          - number of ring structures
        info = {
            "Molecular Formula": Chem.rdMolDescriptors.CalcMolFormula(mol),
            "Molecular Weight": round(float(Descriptors.MolWt(mol)), 3),
            "Heavy Atoms Count": Descriptors.HeavyAtomCount(mol),
            "H Bond Donor Count": Descriptors.NumHDonors(mol),
            "H Bond Acceptor Count": Descriptors.NumHAcceptors(mol),
            "Rotatable Bonds Count": Descriptors.NumRotatableBonds(mol),
            "TPSA": Descriptors.TPSA(mol),
            "Number of Rings": Descriptors.RingCount(mol),
        }
        return info