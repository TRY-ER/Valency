from engine.base import ALLOWED_TYPES

class Visulizer:
    def __init__(self, type: str):
        if type not in ALLOWED_TYPES:
            raise ValueError(f"Invalid visualizer type. type should be one of the {ALLOWED_VISUALIZER_TYPES}")
        self.type = type

    def visualize(self, source_string: str):
        raise NotImplementedError()