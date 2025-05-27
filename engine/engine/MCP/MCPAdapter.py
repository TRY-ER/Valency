import asyncio
from fastmcp.client import Client # Changed from mcp and mcp.client.sse
import json
from mcp.types import TextContent # Assuming TextContent is the correct type for string content

class MCPAdapter:
    def __init__(self, server_url: str):
        """
        Initializes the MCPAdapter with the FastMCP server URL.
        The URL should point to an SSE endpoint. This constructor ensures
        the path ends with "/sse".
        For example:
        - "http://localhost:8080" becomes "http://localhost:8080/sse"
        - "http://localhost:8080/" becomes "http://localhost:8080/sse"
        - "http://localhost:8080/api" becomes "http://localhost:8080/api/sse"
        - "http://localhost:8080/api/" becomes "http://localhost:8080/api/sse"
        - URLs already ending in "/sse" are unchanged.

        Args:
            server_url: The URL of the FastMCP server.
        """
        if server_url.endswith("/sse"):
            self.server_url = server_url
        elif server_url.endswith("/"):
            self.server_url = server_url + "sse"
        else:
            self.server_url = server_url + "/sse"

    async def run_tool(self, tool_name: str, tool_args: dict) -> str:
        """
        Connects to the FastMCP server using fastmcp.client.Client, 
        retrieves a tool by its name, runs the tool with provided arguments,
        and returns the result as a string (often JSON).

        Args:
            tool_name: The name of the tool to run.
            tool_args: A dictionary of arguments for the tool.

        Returns:
            A string containing the tool's result (often JSON), or an error message string.
        
        Raises:
            ConnectionError: If connection to the server fails.
            ValueError: If the tool is not found or if there's an issue with arguments/tool execution.
            Exception: For other underlying errors during the tool call.
        """
        try:
            print(f"Attempting to use FastMCP client with server URL: {self.server_url} for tool '{tool_name}'")

            async with Client(self.server_url) as client:
                # Initialization is handled by Client's context manager
                
                # Validate tool_name
                available_tools_list = await client.list_tools() # Returns list of mcp.common.model.Tool objects
                tool_exists = any(tool.name == tool_name for tool in available_tools_list)

                if not tool_exists:
                    available_tool_names = [tool.name for tool in available_tools_list]
                    error_msg = f"Tool '{tool_name}' not found on the FastMCP server. Available tools: {available_tool_names}"
                    print(error_msg)
                    raise ValueError(error_msg)

                # Validate tool_args type
                if not isinstance(tool_args, dict):
                    # Raising ValueError to be consistent with other argument/tool issues
                    raise ValueError(f"Tool arguments for '{tool_name}' must be a dictionary, got {type(tool_args)}.")

                print(f"Calling tool '{tool_name}' with args: {tool_args} using FastMCP client.")
                result_content = await client.call_tool(tool_name, tool_args)

                print("type of result_content:", type(result_content)) 
                # Ensure the result is a string, serializing dict/list to JSON.
                if isinstance(result_content, list):
                    values = {}
                    for i,d in enumerate(result_content):
                        if isinstance(d, TextContent):
                            # Convert each dict to JSON string
                            print("this is a TextContent")
                            values[i] = d.text
                    return json.dumps(values)
                elif isinstance(result_content, str):
                    return result_content
                else:
                    # For other types like int, float, bool, None convert to string.
                    return str(result_content)

        except ConnectionError as ce: # Covers aiohttp.client_exceptions.ClientConnectorError etc.
            error_msg = f"FastMCP ConnectionError for tool '{tool_name}' at {self.server_url}: {str(ce)}"
            print(error_msg)
            raise ConnectionError(error_msg) # Re-raise standard ConnectionError
        except ValueError as ve: 
            # Catches our "Tool not found" ValueError, arg type ValueError, or ValueErrors from client.call_tool (tool execution error)
            error_msg = f"FastMCP ValueError for tool '{tool_name}': {str(ve)}"
            # Avoid printing if it's our own raised ValueError which already printed.
            if not (str(ve).startswith(f"Tool '{tool_name}' not found") or str(ve).startswith(f"Tool arguments for '{tool_name}' must be a dictionary")):
                 print(error_msg)
            raise # Re-raise the ValueError
        except Exception as e:
            # This could be mcp.common.error.MCPError or other unexpected errors from fastmcp
            error_msg = f"Unexpected FastMCP error for tool '{tool_name}' at {self.server_url}: {type(e).__name__} - {str(e)}"
            print(error_msg)
            raise Exception(error_msg) # Wrap/re-raise as a general Exception with more context
