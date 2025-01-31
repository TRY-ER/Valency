from engine.validator.base import Validator
from rdkit import Chem

class MolValidator(Validator):
    def __init__(self):
        self.type = "MOL"

    def validate(self, source_str: str):
        try:
            mol = Chem.MolFromSmiles(source_str)
            if mol is None:
                return False
            return True
        except:
            return False

if __name__ == "__main__":
    validator = MolValidator()
    print(validator.validate("CC"))
    print(validator.validate("C"))