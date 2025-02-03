class BaseSearcher:
    def query(self, query_val: str, top_k: int = 1):
        raise NotImplementedError("Subclasses must implement this method")