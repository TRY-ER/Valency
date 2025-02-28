import os

def set_env_vars(file_path: str):
    try:
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
    except Exception as e:
        print("Error in loading the .env file. Make sure to add the env files with corresponding paths")