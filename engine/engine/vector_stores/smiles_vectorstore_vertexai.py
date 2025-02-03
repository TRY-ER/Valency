from engine.vector_stores.base import BaseSearcher
from google.cloud import aiplatform_v1
from sentence_transformers import SentenceTransformer

class SMILESVertexAISearcher(BaseSearcher):
    def __init__(self,
                 API_ENDPOINT: str,
                 INDEX_ENDPOINT: str,
                 DEPLOYED_INDEX_ID: str):
        self.API_ENDPOINT = API_ENDPOINT
        self.INDEX_ENDPOINT = INDEX_ENDPOINT
        self.DEPLOYED_INDEX_ID = DEPLOYED_INDEX_ID
        self.model = SentenceTransformer('kuelumbus/polyBERT')

    def get_vector_store_response(self,feature_vector,
                              neighbor_count: int = 1):
        # Configure Vector Search client
        client_options = {"api_endpoint": self.API_ENDPOINT}
        vector_search_client = aiplatform_v1.MatchServiceClient(
            client_options=client_options, )

        # Build FindNeighborsRequest object
        datapoint = aiplatform_v1.IndexDatapoint(feature_vector=feature_vector)

        query = aiplatform_v1.FindNeighborsRequest.Query(
            datapoint=datapoint,

            # The number of nearest neighbors to be retrieved
            neighbor_count=neighbor_count)
        request = aiplatform_v1.FindNeighborsRequest(
            index_endpoint=self.INDEX_ENDPOINT,
            deployed_index_id=self.DEPLOYED_INDEX_ID,
            # Request can have multiple queries
            queries=[query],
            return_full_datapoint=False,
        )

        # Execute the request
        response = vector_search_client.find_neighbors(request)

        # Handle the response
        return response

    def get_embedding(self, smiles: str):
        return self.model.encode(smiles)

    def query(self, query_val: str, top_k: int = 1):
        feature_vector = self.get_embedding(query_val)
        response = self.get_vector_store_response(feature_vector, top_k)
        return response 

if __name__ == "__main__":
    # making the emebddings
    import os

    def set_env_vars(file_path: str):
        with open(file_path, "r") as f:
            lines = f.readlines()
            for line in lines:
                try:
                    line = line.strip()
                    key = line.split("=")[0] 
                    value = line.split("=")[1] 
                    os.environ[key] = value
                except Exception as e:
                    print(f"Error in setting env vars: {e}")
                    continue
    
    set_env_vars("../../../backend/.env")

    SMILES_API_ENDPOINT=os.environ.get["SMILES_API_ENDPOINT"]
    SMILES_INDEX_ENDPOINT=os.environ.get["SMILES_INDEX_ENDPOINT"]
    SMILES_DEPLOYED_INDEX_ID=os.environ.get["SMILES_DEPLOYED_INDEX_ID"]
    PSMILES_API_ENDPOINT=os.environ.get["PSMILES_API_ENDPOINT"]
    PSMILES_INDEX_ENDPOINT=os.environ.get["PSMILES_INDEX_ENDPOINT"]
    PSMILES_DEPLOYED_INDEX_ID=os.environ.get["PSMILES_DEPLOYED_INDEX_ID"]

    smiles_searcher = SMILESVertexAISearcher(PSMILES_API_ENDPOINT, PSMILES_INDEX_ENDPOINT, PSMILES_DEPLOYED_INDEX_ID)
    response = smiles_searcher.query("[*]CCC[*]", 3)
    print("smiles response >>", response)