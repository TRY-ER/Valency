from sentence_transformers import SentenceTransformer
import numpy as np
from typing import List, Dict
import json
from google.cloud import aiplatform
import os

os.environ[
    "GOOGLE_APPLICATION_CREDENTIALS"] = "xx"

PROJECT_ID = "xonext-449315"
LOCATION = "us-central1"
aiplatform.init(project=PROJECT_ID, location=LOCATION)


class TextVectorStore:

    def __init__(self):
        self.model = SentenceTransformer('all-MiniLM-L6-v2')
        self.dimension = 384  # MiniLM dimension

    # def chunk_text(self, text: str, chunk_size: int = 200) -> List[str]:
    #     chunks = text.split("\\n")
    #     return chunks

    def process_file(self, filepath: str, json_path: str):
        with open(filepath, 'r') as f:
            chunks = f.readlines()

        # Create chunks
        # chunks = self.chunk_text(text)

        # Generate embeddings
        embeddings = self.model.encode(chunks)
        data = []
        # Store in Redis
        for idx, (chunk, embedding) in enumerate(zip(chunks, embeddings)):
            embedding = embedding / np.linalg.norm(embedding)
            data.append({
                "id": idx,
                "chunk": chunk.strip(),
                "embedding": embedding.tolist()
            })

        # Write to JSON
        self.write_to_json(data, json_path)

        return len(chunks)

    def write_to_json(self, data: List[Dict], filename: str):
        """Write data to a JSON file in a single line"""
        with open(filename, 'w') as f:
            for record in data:
                line = json.dumps(record, separators=(',', ':'))
                f.write(line + '\n')


def get_vector_store_response(feature_vector,
                              neighbor_count: int = 1):
    from google.cloud import aiplatform_v1

    # Set variables for the current deployed index.
    API_ENDPOINT = "xx"
    INDEX_ENDPOINT = "xx"
    DEPLOYED_INDEX_ID = "xx"

    # Configure Vector Search client
    client_options = {"api_endpoint": API_ENDPOINT}
    vector_search_client = aiplatform_v1.MatchServiceClient(
        client_options=client_options, )

    # Build FindNeighborsRequest object
    datapoint = aiplatform_v1.IndexDatapoint(feature_vector=feature_vector)

    query = aiplatform_v1.FindNeighborsRequest.Query(
        datapoint=datapoint,

        # The number of nearest neighbors to be retrieved
        neighbor_count=neighbor_count)
    request = aiplatform_v1.FindNeighborsRequest(
        index_endpoint=INDEX_ENDPOINT,
        deployed_index_id=DEPLOYED_INDEX_ID,
        # Request can have multiple queries
        queries=[query],
        return_full_datapoint=True,
    )

    # Execute the request
    response = vector_search_client.find_neighbors(request)

    # Handle the response
    return response


if __name__ == "__main__":
    # making the emebddings

    # store = TextVectorStore()
    # filepath = "sample_vector_store_content.txt"  # Create this file with sample text
    # json_path = "./sample_vector_store_content.json"
    # chunk_count = store.process_file(filepath, json_path=json_path)

    # setting up the cloud vectro store with index
    # bucket_uri = "gs://sample-vector-bucket"

    # index = aiplatform.MatchingEngineIndex.create_tree_ah_index(
    #     display_name="sample-vector-index",
    #     contents_delta_uri=bucket_uri,
    #     dimensions=384,
    #     approximate_neighbors_count=10
    # )

    # print("index >>", index)
    query_str = "how to prepare a pizza ?"
    model = SentenceTransformer('all-MiniLM-L6-v2')
    query_vector = model.encode(query_str) 
    response = get_vector_store_response(feature_vector=query_vector)
    print("response >>", response)
