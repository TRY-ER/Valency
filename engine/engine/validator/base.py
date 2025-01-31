from engine.base import ALLOWED_TYPES

class Validator:
    def __init__(type: str):
        if type not in ALLOWED_TYPES:
            raise ValueError(f"Type {type} is not allowed, Try one of {ALLOWED_TYPES}")
        self.type = type 

    def validate(source_str: str):
        raise NotImplementedError("This method is not implemented. It should be implemented in child classes")