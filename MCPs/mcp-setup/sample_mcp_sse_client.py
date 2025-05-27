# sse_client.py
import asyncio
from fastmcp.client import Client # Corrected import based on jlowin/fastmcp documentation

async def main():
    # Server URL for FastMCP (base URL, SSE endpoint is discovered)
    # The sample server was set to run on port 8080.
    server_url = "http://127.0.0.1:8008/sse"
    
    print(f"Connecting to FastMCP server at {server_url}...")
    
    # Create the Client session
    async with Client(server_url) as client: # Changed to Client from fastmcp.client
        # Initialization is handled by Client context manager
        
        # List available tools
        tools = await client.list_tools()
        # The tool objects themselves might not directly have a .name attribute in a simple list.
        # fastmcp client's list_tools() returns a list of mcp.common.model.Tool objects.
        # These objects do have a .name attribute.
        print("Available tools:", [tool.name for tool in tools])
        
        # Call the greet tool
        # The response from call_tool in fastmcp is the direct result.
        greeting_result = await client.call_tool("greet", {"name": "Bob"})
        print("Greeting result:", greeting_result)
        
        # Call the add tool
        addition_result = await client.call_tool("add", {"a": 10, "b": 32})
        print("Addition result:", addition_result)

if __name__ == "__main__":
    asyncio.run(main())