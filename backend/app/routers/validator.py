from fastapi import APIRouter
from engine.validator.validators import MolValidator, ProtValidator, PolymerValidator
from engine.base import ALLOWED_TYPES

router = APIRouter()

validator_mapper = {
    "MOL": MolValidator(),
    "PROT": ProtValidator(),
    "POLY": PolymerValidator()
}

@router.post("/")
def validate_molecule(data: dict):
    validate_type = data["type"]
    if validate_type is None or validate_type == "" or validate_type not in ALLOWED_TYPES:
        return {"status": "error",
                "message": f"Validation type should be one of {ALLOWED_TYPES}"}
    validator = validator_mapper.get(validate_type)
    return {"status": "success",
            "valid": validator.validate(str(data["value"]))}

