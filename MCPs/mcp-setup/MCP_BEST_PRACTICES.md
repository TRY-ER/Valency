# MCP Server Development Best Practices

This document outlines standard practices to follow when developing MCP (Model Context Protocol) servers to ensure consistency, stability, and ease of configuration.

## 1. Standardized Tool Output: JSON Strings

All tools defined within an MCP server should consistently return their output as a **JSON string**. This applies to both successful results and error conditions.

### Rationale:
- **Uniformity:** Provides a predictable output format for any client or tool consuming the MCP.
- **Robustness:** JSON is a widely supported, language-agnostic data interchange format.
- **Error Handling:** Structuring errors as JSON objects within the string allows for detailed and parsable error information.

### Implementation:

#### a. Successful Data:
- Wrap the primary data (e.g., a list of results, a dictionary of details) within a main key in a dictionary (e.g., `"data": ...`, `"count": ...`, or a more specific key like `"svg_image": ...`).
- Use `json.dumps()` to serialize this dictionary into a JSON string.
- Update the tool function's return type hint to `str`.

**Example (Python):**
```python
import json
from mcp.server.fastmcp import FastMCP # Assuming FastMCP is used

mcp = FastMCP("MyServer")

@mcp.tool()
def get_some_data(param: str) -> str: # Return type is str
    try:
        # ... your logic to fetch or compute data ...
        actual_data = [{"id": 1, "value": "example"}, {"id": 2, "value": "another"}]
        
        # Wrap data and serialize to JSON string
        return json.dumps({"data": actual_data}) 
    except Exception as e:
        # ... (see error handling below) ...
        error_response = {"error": "Failed to process request", "details": str(e)}
        return json.dumps(error_response)

```

#### b. Error Conditions:
- Catch exceptions within your tool's logic.
- Construct a dictionary containing at least an `"error"` key with a descriptive message. Optionally, include a `"details"` key for more specific error information (like the exception string).
- Use `json.dumps()` to serialize this error dictionary into a JSON string.

**Example (Python):**
```python
# ... (inside a tool function) ...
try:
    # ... risky operation ...
    if some_condition_fails:
        raise ValueError("A specific condition was not met.")
    # ...
    return json.dumps({"data": "Success!"})
except ValueError as ve:
    return json.dumps({"error": "Invalid input provided.", "details": str(ve)})
except requests.exceptions.HTTPError as http_e: # Example for network errors
    err_text = http_e.response.text if http_e.response else "No response text"
    return json.dumps({"error": f"HTTP error: {http_e}", "details": err_text})
except Exception as e: # Generic fallback
    return json.dumps({"error": "An unexpected error occurred.", "details": str(e)})
```

## 2. Configuration via `.env` Files

For managing server configurations like host and port, use `.env` files. This practice keeps sensitive or environment-specific settings out of version control and allows for easy modification without code changes.

### Rationale:
- **Security:** Prevents hardcoding sensitive information.
- **Flexibility:** Allows different configurations for development, staging, and production environments.
- **Standardization:** `.env` files are a common convention.

### Implementation:

#### a. Create/Update `.env` File:
In the root directory of your MCP server project (or a designated config location), create or update a `.env` file. For each MCP server, define its host and port. Use a consistent naming convention, for example:
```env
# For RCSB PDB MCP Server
RCSB_HOST=0.0.0.0
RCSB_PORT=8052

# For ChEMBL MCP Server
CHEMBL_HOST=0.0.0.0
CHEMBL_PORT=8051

# For Deep-PK MCP Server
DEEP_PK_MCP_HOST=0.0.0.0
DEEP_PK_MCP_PORT=8053

# For general MCP settings if applicable
# MCP_TRANSPORT=sse # or stdio
```

#### b. Load Environment Variables in Your Server Script:
Use a library like `python-dotenv` (or your custom `env_loader.py`) at the beginning of your server script to load these variables from the `.env` file into the environment.

**Example (Python using `python-dotenv`):**
```python
import os
import json
from dotenv import load_dotenv
from mcp.server.fastmcp import FastMCP

# Load environment variables from .env file
load_dotenv() 

# Configuration for this specific MCP Server (e.g., Deep-PK)
DEFAULT_MCP_HOST = "0.0.0.0" # Default if not in .env
DEFAULT_MCP_PORT = 8053   # Default if not in .env

# Get host and port from environment variables, with fallbacks
mcp_host = os.getenv("DEEP_PK_MCP_HOST", DEFAULT_MCP_HOST)
mcp_port = int(os.getenv("DEEP_PK_MCP_PORT", str(DEFAULT_MCP_PORT)))

mcp = FastMCP(
    "My MCP Server Name",
    description="Description of my server.",
    host=mcp_host,
    port=mcp_port
)

# ... rest of your server code and tool definitions ...

if __name__ == "__main__":
    print(f"Starting My MCP Server on {mcp_host}:{mcp_port}")
    mcp.run() # Add transport="sse" or "stdio" as needed
```

**Example (Python using a custom `env_loader.py`):**
```python
# In your server_script.py
import os
import json
from env_loader import load_env_vars # Assuming your custom loader
from mcp.server.fastmcp import FastMCP

load_env_vars() # Call your custom loader

# ... rest of the configuration is similar to the python-dotenv example ...
```

#### c. Add `.env` to `.gitignore`:
Ensure that the `.env` file itself is **not** committed to your version control system (e.g., Git). Add it to your `.gitignore` file:
```
.env
```
Provide a template or example file (e.g., `.env.example`) in your repository that developers can copy to create their local `.env` file.

## Summary
By adhering to these practices, MCP server development becomes more standardized, robust, and easier to manage across different environments and by multiple developers.
