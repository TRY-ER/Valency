
from google.adk.tools.base_tool import BaseTool
from google.adk.tools.tool_context import ToolContext
from typing import Optional, Dict, Any, List
from mcp.types import CallToolResult, TextContent
import json
import uuid

from custom_utils.mongodb_handler import MongoDBHandler

MAX_LENGTH = 2000  # Define a maximum length for tool output

def after_tool_output_limit_callback(tool: BaseTool, args: Dict[str, Any], tool_context: ToolContext, tool_response: CallToolResult):
    """
    Callback function to handle tool output overflow.
    
    This function checks if the tool response exceeds a certain length and truncates it if necessary.
    It also stores the complete tool response in MongoDB for later retrieval.
    
    Args:
        tool (BaseTool): The tool that produced the response.
        args (Dict[str, Any]): The arguments passed to the tool.
        tool_context (ToolContext): The context in which the tool was executed.
        tool_response (CallToolResult): The response from the tool.
        
    Returns:
        Optional[dict]: The modified tool response, or None if no modification is needed.
    """
    # print("tool_response", tool_response)
    # print("tool_response type", type(tool_response))
    
    # Extract session information from tool_context
    user_id = getattr(tool_context, 'user_id', 'unknown_user')
    session_id = getattr(tool_context, 'session_id', 'unknown_session')
    
    # Generate a unique tool_id
    db_handler = MongoDBHandler()
    tool_id = db_handler.generate_tool_id() 
    
    # Process the tool response based on its type
    new_tool_responses = []
    response_is_truncated = False
    
    if isinstance(tool_response, CallToolResult):
        print("tool_response is CallToolResult")
        
        # Store the complete response in MongoDB before truncating
        try:
            # Convert the tool response to a serializable format
            serializable_response = {
                "isError": tool_response.isError,
                "content": [{"type": c.type, "text": c.text} for c in tool_response.content] if hasattr(tool_response, 'content') else []
            }
            
            # Store in MongoDB
            mongo_id = db_handler.store_tool_response(
                tool_id=tool_id,
                tool_name=tool.__class__.__name__,
                user_id=user_id,
                session_id=session_id,
                response_data=serializable_response
            )
            print(f"Stored complete tool response in MongoDB with ID: {mongo_id}")
        except Exception as e:
            print(f"Error storing tool response in MongoDB: {str(e)}")
        
        # Process and truncate if needed
        if not tool_response.isError:
            print("tool_response is not error")
            content = tool_response.content
            for c in content:
                text_content = c.text
                if len(text_content) > MAX_LENGTH:
                    truncated_text = text_content[:MAX_LENGTH] + f'... [truncated]'
                    response_is_truncated = True
                else:
                    truncated_text = text_content
                
                new_text_response = TextContent(type="text", text=truncated_text)
                new_tool_responses.append(new_text_response)
            
            # If response was truncated, add a special message with the tool_id
            # if response_is_truncated:
            #     retrieval_info = TextContent(
            #         type="text", 
            #         text=f"Note: Some tool output was truncated. Full response can be retrieved using tool_id: {tool_id}"
            #     )
            #     new_tool_responses.append(retrieval_info)
            
            new_tool_result = CallToolResult(
                content=new_tool_responses,
                isError=False
            )
            new_tool_result.meta = {'id': tool_id, 'tool_name': tool.__class__.__name__}
            return new_tool_result
        else:
            return tool_response
    elif isinstance(tool_response, str):
        # If the response is a string, store it in MongoDB and truncate if necessary
        try:
            # Store the string response in MongoDB
            mongo_id = db_handler.store_tool_response(
                tool_id=tool_id,
                tool_name=tool.__class__.__name__,
                user_id=user_id,
                session_id=session_id,
                response_data={"text": tool_response}
            )
            print(f"Stored string tool response in MongoDB with ID: {mongo_id}")
        except Exception as e:
            print(f"Error storing string tool response in MongoDB: {str(e)}")
            
        # Truncate if needed
        if len(tool_response) > MAX_LENGTH:
            return tool_response[:MAX_LENGTH] + f'... [truncated, full response available with tool_id: {tool_id}]'
        else:
            return tool_response
    else:
        # For other types, attempt to store in MongoDB but return unchanged
        try:
            mongo_id = db_handler.store_tool_response(
                tool_id=tool_id,
                tool_name=tool.__class__.__name__,
                user_id=user_id,
                session_id=session_id,
                response_data=tool_response
            )
            print(f"Stored other type tool response in MongoDB with ID: {mongo_id}")
        except Exception as e:
            print(f"Error storing other type tool response in MongoDB: {str(e)}")
            
        return tool_response