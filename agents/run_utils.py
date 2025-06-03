from google.genai import types # Added for types.FunctionResponse
# from google.adk import CallToolResult
import json
from mcp.types import CallToolResult

# ANSI color codes for terminal output
class Colors:
    RESET = "\033[0m"
    BOLD = "\033[1m"
    UNDERLINE = "\033[4m"

    # Foreground colors
    BLACK = "\033[30m"
    RED = "\033[31m"
    GREEN = "\033[32m"
    YELLOW = "\033[33m"
    BLUE = "\033[34m"
    MAGENTA = "\033[35m"
    CYAN = "\033[36m"
    WHITE = "\033[37m"

    # Background colors
    BG_BLACK = "\033[40m"
    BG_RED = "\033[41m"
    BG_GREEN = "\033[42m"
    BG_YELLOW = "\033[43m"
    BG_BLUE = "\033[44m"
    BG_MAGENTA = "\033[45m"
    BG_CYAN = "\033[46m"
    BG_WHITE = "\033[47m"


async def display_state(
    session_service, app_name, user_id, session_id, label="Current State"
):
    """Display the current session state in a formatted way."""
    try:
        session = await session_service.get_session(
            app_name=app_name, user_id=user_id, session_id=session_id
        )

        # Format the output with clear sections
        print(f"\n{'-' * 10} {label} {'-' * 10}")

        # Handle the user name
        user_name = session.state.get("user_name", "Unknown")
        print(f"👤 User: {user_name}")

        # Handle reminders
        reminders = session.state.get("reminders", [])
        if reminders:
            print("📝 Reminders:")
            for idx, reminder in enumerate(reminders, 1):
                print(f"  {idx}. {reminder}")
        else:
            print("📝 Reminders: None")

        print("-" * (22 + len(label)))
    except Exception as e:
        print(f"Error displaying state: {e}")


async def process_agent_response(event):
    """Process agent response events and yield structured Python dictionaries or specific strings."""
    # Log basic event info
    print(f"Event ID: {event.id}, Author: {event.author}")
    function_response_id = None

    if event.content and event.content.parts:
        for part in event.content.parts:
            if hasattr(part, "function_call") and part.function_call:
                function_data = {
                    "name": part.function_call.name,
                    "id": part.function_call.id,
                    # Ensure args are serializable if they exist
                    "args": part.function_call.args if hasattr(part.function_call, "args") else {}
                }
                data = { "type": "function_call", "data": function_data}
                yield data
            
            elif hasattr(part, "function_response") and part.function_response:
                # Print tool response information
                tool_response = {} 
                if 'result' in part.function_response.response:
                    result = part.function_response.response['result']
                    # Ensure result is serializable; ADK's FunctionResponse might need conversion
                    print("result >>", result)
                    print("result type >>", type(result))
                    if isinstance(result, CallToolResult): # Assuming types.FunctionResponse is from google.genai
                        if result.isError:
                            tool_response = {
                                "type": "error",
                                "text_content": "" # Or some error message from result if available
                            }
                        else:
                            tool_contents = ""
                            # Assuming result.content is iterable and parts have .text
                            for r_part in result.content: # Renamed 'r' to 'r_part' to avoid conflict if 'r' is used above
                                if hasattr(r_part, 'text'): # Check if it's a text part
                                    tool_contents += r_part.text
                                # Handle other part types if necessary, e.g. images, etc.
                            
                            tool_response = {
                                "type": "success",
                                "text_content": tool_contents
                            }
                        if result.meta and 'id' in result.meta:
                            function_response_id = result.meta['id']
                    else: # Handle cases where result is not a FunctionResponse object (e.g., already a dict or string)
                        tool_response = {
                            "type": "unknown", # Or try to infer type
                            "text_content": str(result) # Fallback to string representation
                        }

                response_data = {
                    "name": part.function_response.name,
                    "id": part.function_response.id,
                    "response": tool_response,
                    "function_response_id": function_response_id 
                }
                data = { "type": "function_response", "data": response_data}
                yield data

            # Also print any text parts found in any event for debugging
            elif hasattr(part, "text") and part.text and not part.text.isspace():
                # print(f"  Text: '{part.text.strip()}'")
                yield {"type": "text", "content": part.text.strip()}

    # Check for final response after specific parts
    final_response = None
    if event.is_final_response():
        if (
            event.content
            and event.content.parts
            and hasattr(event.content.parts[0], "text")
            and event.content.parts[0].text
        ):
            final_response = event.content.parts[0].text.strip()
            print(
                f"\n{Colors.BG_BLUE}{Colors.WHITE}{Colors.BOLD}╔══ AGENT RESPONSE ═════════════════════════════════════════{Colors.RESET}"
            )
            print(f"{Colors.CYAN}{Colors.BOLD}{final_response}{Colors.RESET}")
            print(
                f"{Colors.BG_BLUE}{Colors.WHITE}{Colors.BOLD}╚═════════════════════════════════════════════════════════════{Colors.RESET}\n"
            )
            # yield {"type": "text", "content": final_response}
            yield "<|done|>" # Signal completion
        else:
            print(
                f"\n{Colors.BG_RED}{Colors.WHITE}{Colors.BOLD}==> Final Agent Response: [No text content in final event]{Colors.RESET}\n"
            )
            # Still yield done if it's a final response without text, to close stream.
            yield "<|done|>"

async def call_agent_async(runner, user_id, session_id, query):
    """Call the agent asynchronously with the user's query."""
    content = types.Content(role="user", parts=[types.Part(text=query)])
    print(
        f"\n{Colors.BG_GREEN}{Colors.BLACK}{Colors.BOLD}--- Running Query: {query} ---{Colors.RESET}"
    )
    final_response_text = None

    # # Display state before processing
    # await display_state(
    #     runner.session_service,
    #     runner.app_name,
    #     user_id,
    #     session_id,
    #     "State BEFORE processing",
    # )

    try:
        async for event in runner.run_async(
            user_id=user_id, session_id=session_id, new_message=content
        ):
            # Process each event and get the final response if available
            async for event_response in process_agent_response(event):
                yield event_response 
            # if response:
            #     final_response_text = response
    except Exception as e:
        print(f"Error during agent call: {e}")

    # # Display state after processing the message
    # await display_state(
    #     runner.session_service,
    #     runner.app_name,
    #     user_id,
    #     session_id,
    #     "State AFTER processing",
    # )

    # return final_response_text

def parse_event_list(event_list: list) -> list:
    """
    Schema required at the frontend
            [{
                "query": <query_text>,
                "response": [<resopnse_events_as_is>, ...]
            }] 
    """
    parsed_events = []
    user_query = None
    for i, event in enumerate(event_list):
        if event.content and event.content.parts:
            if event.content.role == "user":
                if event.content.parts[0].text and not event.content.parts[0].text.isspace():
                    if user_query is not None and "response" in user_query:
                        if len(user_query["response"]) > 0:
                            # print("user query added !")
                            parsed_events.append(user_query)
                    user_query = {"query": event.content.parts[0].text.strip(), "response": []}
                for part in event.content.parts:
                    if hasattr(part, "function_response") and part.function_response:
                            function_response_id = None
                            tool_response = {} 
                            if 'result' in part.function_response.response:
                                result = part.function_response.response['result']
                                if isinstance(result, dict):
                                    if result["isError"]:
                                        tool_response = {
                                            "type": "error",
                                            "text_content": "" # Or some error message from result if available
                                        }
                                    else:
                                        tool_contents = ""
                                        for r_part in result["content"]:
                                            if hasattr(r_part, 'text'):
                                                tool_contents += r_part.text
                                        tool_response = {
                                            "type": "success",
                                            "text_content": tool_contents
                                        }
                                    if result["meta"] and 'id' in result["meta"]:
                                        function_response_id = result["meta"]['id']
                                else:
                                    tool_response = {
                                        "type": "unknown",
                                        "text_content": str(result) # Fallback to string representation
                                    }
                            response_data = {
                                "name": part.function_response.name,
                                "id": part.function_response.id,
                                "response": tool_response,
                                "function_response_id": function_response_id
                            }
                            if user_query is not None:
                                user_query["response"].append({"type": "function_response", "data": response_data})
            if event.content.role == "model":
                if user_query is not None:
                    # Process the event parts and yield structured data
                    for part in event.content.parts:
                        if hasattr(part, "text") and part.text and not part.text.isspace():
                            # print("text condition in parse event in model response >>", part.text)
                            user_query["response"].append({"type": "text", "content": part.text.strip()})
                        elif hasattr(part, "function_call") and part.function_call:
                            function_data = {
                                "name": part.function_call.name,
                                "id": part.function_call.id,
                                "args": part.function_call.args if hasattr(part.function_call, "args") else {}
                            }
                            user_query["response"].append({"type": "function_call", "data": function_data})
    if user_query is not None and "response" in user_query:
        parsed_events.append(user_query)
    # print('parsed_events >>', parsed_events)
    return parsed_events 