import json
import chromadb

def create_chroma_store_from_json(json_path: str, collection_name: str = "smiles_data", persist_directory: str = "./chroma_store"):
    # Initialize Chroma
    client = chromadb.PersistentClient(
        path=persist_directory
    )
    collection = client.get_or_create_collection(name=collection_name)

    # Read the JSON file (assuming NDJSON with one object per line)
    with open(json_path, 'r') as f:
        for line in f:
            if not line.strip():
                continue
            record = json.loads(line)
            doc_id = str(record["id"])
            embedding = record["embedding"]
            metadata = record["data"]

            collection.add(
                documents=[metadata.get("smiles", doc_id)],
                embeddings=[embedding],
                metadatas=[metadata],
                ids=[doc_id]
            )

    print("Chroma store created successfully!")

if __name__ == "__main__":
    # setting up vector store for the psmiles data:

    json_path = "../resources/chroma_store/input_data/psmiles.json"
    collection_name = "psmiles_data"
    chroma_store_path = "../resources/chroma_store/vector_stores" 
    create_chroma_store_from_json(json_path, collection_name, persist_directory=chroma_store_path)

    # setting up vector store for the smiles data:

    json_path = "../resources/chroma_store/input_data/smiles.json"
    collection_name = "smiles_data"
    chroma_store_path = "../resources/chroma_store/vector_stores" 
    create_chroma_store_from_json(json_path, collection_name, persist_directory=chroma_store_path)

    # Query a sample embedding

    client = chromadb.PersistentClient(
        path=chroma_store_path
    )
    collection = client.get_collection(name=collection_name)

    # # A dummy embedding vector for demonstration (adjust to match the actual model dimension)
    sample_embedding = [0.0] * 600
    result = collection.query(query_embeddings=[sample_embedding], n_results=1)
    print("Query result:", result)

    # querying a sample psmiles

    # collection_name = "psmiles_data"
    # chroma_store_path = "../../../bulk_data/chroma_store" 
 
    # from sentence_transformers import SentenceTransformer
    # model = SentenceTransformer('kuelumbus/polyBERT')
    # client = chromadb.PersistentClient(
    #     path=chroma_store_path
    # )
    # q_vec = model.encode("[*]CCO[*]") 
    # collection = client.get_collection(name=collection_name)
    # result = collection.query(query_embeddings=[q_vec], n_results=3)
    # print("Query result:", result)


    # querying a sample smiles

    # collection_name = "smiles_data"
    # chroma_store_path = "../../../bulk_data/chroma_store" 
 
    # from sentence_transformers import SentenceTransformer
    # model = SentenceTransformer('kuelumbus/polyBERT')
    # client = chromadb.PersistentClient(
    #     path=chroma_store_path
    # )
    # q_vec = model.encode("CCO") 
    # collection = client.get_collection(name=collection_name)
    # result = collection.query(query_embeddings=[q_vec], n_results=3)
    # print("Query result:", result)