import chromadb
from chromadb.config import Settings
from sentence_transformers import SentenceTransformer
from engine.vector_stores.base import BaseSearcher


class ChromaSearcher(BaseSearcher):
    def __init__(self,
                 collection_name: str,
                 persist_directory: str):
        self.persist_directory = persist_directory
        self.client = chromadb.PersistentClient(path=self.persist_directory)
        self.collection = self.client.get_or_create_collection(
            name=collection_name)
        self.model = SentenceTransformer('kuelumbus/polyBERT')

    def get_embedding(self, smiles: str):
        return self.model.encode(smiles).tolist()

    def query(self, query_val: str, top_k: int = 3):
        embedding = self.get_embedding(query_val)
        return self.collection.query(query_embeddings=[embedding],
                                     n_results=top_k)