import os

def load_env_vars(env_path=".env"):
    """Loads environment variables from a .env file line by line."""
    try:
        with open(env_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line and '=' in line and not line.startswith('#'):
                    key, value = line.split('=', 1)
                    key = key.strip()
                    value = value.strip()
                    # Remove surrounding quotes if any
                    if (value.startswith(''') and value.endswith(''')) or \
                       (value.startswith('"') and value.endswith('"')):
                        value = value[1:-1]
                    os.environ[key] = value
    except FileNotFoundError:
        print(f"Warning: {env_path} file not found. Skipping environment variable loading from file.")
    except Exception as e:
        print(f"Error loading {env_path}: {e}")
