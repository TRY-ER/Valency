
from google.adk.tools.base_tool import BaseTool
from google.adk.tools.tool_context import ToolContext
from typing import Optional, Dict, Any
from mcp.types import CallToolResult, TextContent


MAX_LENGTH = 2000  # Define a maximum length for tool output

def after_tool_output_limit_callback(tool: BaseTool, args: Dict[str, Any], tool_context: ToolContext, tool_response: CallToolResult) -> Optional[dict]:
    """
    Callback function to handle tool output overflow.
    
    This function checks if the tool response exceeds a certain length and truncates it if necessary.
    
    Args:
        tool (BaseTool): The tool that produced the response.
        tool_context (ToolContext): The context in which the tool was executed.
        tool_response (dict): The response from the tool.
        
    Returns:
        Optional[dict]: The modified tool response, or None if no modification is needed.
    """
    print("tool_response", tool_response)
    print("tool_response type", type(tool_response))
    new_tool_responses = []
    if isinstance(tool_response, CallToolResult):
        print("tool_response is CallToolResult") # Debugging line
        if not tool_response.isError:
            print("tool_response is not error") # Debugging line
            content = tool_response.content
            for c in content:
                tool_response = c.text
                tool_response = tool_response[:MAX_LENGTH] + '... [truncated]' if len(tool_response) > MAX_LENGTH else tool_response
                new_text_response = TextContent(type="text", text=tool_response)
                new_tool_responses.append(new_text_response) 
            new_tool_result = CallToolResult(
                content=new_tool_responses,
                isError=False)
            return new_tool_result
        else:
            return tool_response
    elif isinstance(tool_response, str):
        # If the response is a string, truncate it if necessary
        tool_response = tool_response[:MAX_LENGTH] + '... [truncated]' if len(tool_response) > MAX_LENGTH else tool_response
        return tool_response
    else:
        return tool_response