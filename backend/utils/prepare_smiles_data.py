from sentence_transformers import SentenceTransformer
import pandas as pd
import json

# encoded = polyBERT.encode("CCC")

class SMILESVectorGenrator_VertexAI:
    def __init__(self):
        self.model = SentenceTransformer('kuelumbus/polyBERT')
        self.dimension = 600

    def process_csv(self, filepath: str, json_path: str, verbose_interval = 100):
        df = pd.read_csv(filepath)
        temp = set()
        with open(filepath, 'r') as f:
            for index, row in df[:10].iterrows():
                source_string = row['SMILES']
                if source_string in temp:
                    continue
                temp.add(source_string)
                data = {}
                for i,j in row.items():
                    if i != "mol":
                        data[i] = j
                embeddings = self.model.encode(source_string) 
                payload = {
                    "id": index,
                    "data": data,
                    "embedding": embeddings.tolist()
                }
                line = json.dumps(payload, separators=(',', ':'))
                with open(json_path, 'a') as f:
                    f.write(line + '\n')
                if index % verbose_interval == 0:
                    print(f"[++] Processed {index} records")

class PSMILESVectorGenrator_VertexAI:
    def __init__(self):
        self.model = SentenceTransformer('kuelumbus/polyBERT')
        self.dimension = 600

    def process_csv(self, filepath: str, json_path: str, verbose_interval = 100):
        df = pd.read_csv(filepath)
        temp = set()
        with open(filepath, 'r') as f:
            for index, row in df[:10].iterrows():
                source_string = row['smiles']
                if source_string in temp:
                    continue
                temp.add(source_string)
                data = {}
                for i,j in row.items():
                    if i != "mol":
                        data[i] = j
                embeddings = self.model.encode(source_string) 
                payload = {
                    "id": index,
                    "data": data,
                    "embedding": embeddings.tolist()
                }
                line = json.dumps(payload, separators=(',', ':'))
                with open(json_path, 'a') as f:
                    f.write(line + '\n')
                if index % verbose_interval == 0:
                    print(f"[++] Processed {index} records")

if __name__ == "__main__":
    # vg = SMILESVectorGenrator()
    # vg.process_csv("data/SMILES_Big_Data_Set.csv", "output/smiles.json")
    vg = PSMILESVectorGenrator_VertexAI()
    vg.process_csv("data/psmiles_source/bandgap_chain.csv", "output/psmiles.json")