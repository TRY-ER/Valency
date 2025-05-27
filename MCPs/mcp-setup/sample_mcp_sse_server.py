# sse_server.py
from fastmcp import FastMCP # Changed import

# Create the MCP server
# Default port for fastmcp is 8000, but we'll set it to 8080 to match the original example for now.
mcp = FastMCP(
    "SSE Example Server",
    settings={"host": "0.0.0.0", "port": "3003"} 
)

@mcp.tool()
def greet(name: str) -> str:
    """Greet a user by name"""
    return f"Hello, {name}! Welcome to the SSE server."

@mcp.tool()
def add(a: int, b: int) -> str:
    """Add two numbers and return the result"""
    return f"The sum of {a} and {b} is {a + b}."

if __name__ == "__main__":
    port = mcp.settings.port # Get port from mcp settings
    host = mcp.settings.host # Get host from mcp settings

    print(f"Starting MCP server with FastMCP on port {port}...")
    # FastMCP typically exposes SSE at /mcp/sse by default
    # The exact URL might depend on fastmcp's internal routing,
    # but the client will connect to the base URL and discover endpoints.
    print(f"MCP server (FastAPI) running at: http://{host}:{port}")
    print(f"Default SSE endpoint likely at: http://{host}:{port}/mcp/sse") # Informational
        
    mcp.run(transport="sse", host="127.0.0.1", port=8008) # fastmcp handles uvicorn internally