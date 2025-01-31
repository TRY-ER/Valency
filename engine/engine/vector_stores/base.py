from abc import ABC, abstractmethod
import numpy as np
from typing import List, Dict, Any, Optional
import redis
from redis.commands.search.field import VectorField
from redis.commands.search.query import Query

class VectorStore(ABC):
    def __init__(self, dimension: int, redis_host: str = 'localhost', redis_port: int = 6379):
        """Initialize vector store with given dimension"""
        self.dimension = dimension
        self.redis_client = redis.Redis(host=redis_host, port=redis_port)
        self.index_name = "vector_idx"
        
        # Create Redis search index if it doesn't exist
        try:
            self.redis_client.ft(self.index_name).info()
        except:
            # Define the vector field
            vector_field = VectorField("vector",
                                     "FLAT",
                                     {
                                         "TYPE": "FLOAT32",
                                         "DIM": self.dimension,
                                         "DISTANCE_METRIC": "COSINE"
                                     })
            # Create the index
            self.redis_client.ft(self.index_name).create_index([vector_field])
    
    def add_vectors(self, vectors: List[np.ndarray], metadata: Optional[List[Dict[str, Any]]] = None):
        """Add vectors with optional metadata to store"""
        if metadata is None:
            metadata = [{} for _ in vectors]
            
        for idx, (vector, meta) in enumerate(zip(vectors, metadata)):
            # Normalize vector
            vector = vector / np.linalg.norm(vector)
            
            # Prepare vector data
            key = f"vector:{idx}"
            vector_data = {
                "vector": vector.astype(np.float32).tobytes(),
                **meta
            }
            
            # Store in Redis
            self.redis_client.hset(key, mapping=vector_data)
    
    def search_similar(self, query_vector: np.ndarray, k: int = 5) -> List[Dict[str, Any]]:
        """Search k most similar vectors"""
        # Normalize query vector
        query_vector = query_vector / np.linalg.norm(query_vector)
        
        # Prepare search query
        q = Query(f"*=>[KNN {k} @vector $query_vector AS score]")\
            .dialect(2)\
            .return_fields("score", "__key")\
            .paging(0, k)\
            .sort_by("score")
            
        # Execute search
        query_params = {"query_vector": query_vector.astype(np.float32).tobytes()}
        results = self.redis_client.ft(self.index_name).search(q, query_params)
        
        # Format results
        similar_vectors = []
        for doc in results.docs:
            similar_vectors.append({
                "id": doc.id,
                "score": float(doc.score),
                "metadata": self.redis_client.hgetall(doc.id)
            })
            
        return similar_vectors

    def __len__(self):
        """Return number of vectors in store"""
        return self.redis_client.dbsize()

if __name__ == "__main__":
    # Example usage
    vector_store = VectorStore(dimension=768)  # For BERT embeddings

    # Add vectors
    vectors = [np.random.randn(768) for _ in range(10)]
    metadata = [{"text": f"Document {i}"} for i in range(10)]
    vector_store.add_vectors(vectors, metadata)

    # Search similar
    query = np.random.randn(768)
    results = vector_store.search_similar(query, k=3)