from rdkit import Chem
from rdkit.Chem import Descriptors

class PolymerInfo:
    def get_info(self, source_str: str) -> dict:
        """
        Return basic RDKit descriptors from a SMILES string.
        """
        mol = Chem.MolFromSmiles(source_str)

        # get indexes of [*] content in the PSMILES string
        wildcard_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() == '*']
        wildcard_indices = ",".join([str(index) for index in wildcard_indices])

        info = {
            "Molecular Formula": Chem.rdMolDescriptors.CalcMolFormula(mol),
            "Monomer Molecular Weight": round(float(Descriptors.MolWt(mol)), 3),
            "Number of Rings in Monomer": Descriptors.RingCount(mol),
            "Open Bond Indexes": wildcard_indices
        }
        
        return info


if __name__ == "__main__":
    informer = PolymerInfo()
    print(informer.get_info("[*]CC[*]"))  # Valid PSMILES