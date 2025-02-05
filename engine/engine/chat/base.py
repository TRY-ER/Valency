class BaseLLM:
    def respond(self, query, *kwargs):
        raise NotImplementedError("This method should be implemented by the subclass")

    def stream_respond(self, query, *kwargs):
        raise NotImplementedError("This method should be implemented by the subclass")

class BaseFormater:
    def format(self, response, *kwargs):
        raise NotImplementedError("This method should be implemented by the subclass")