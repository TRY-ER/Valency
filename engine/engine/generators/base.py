class Generator:
    def __init__(self, input_types: list):
        self.input_types = input_type 

    def generate(self):
        raise NotImplementedError("This method should be implemented in the subclass")

    def stream_generate(self):
        raise NotImplementedError("This method should be implemented in the subclass")