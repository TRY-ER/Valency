from typing import List, Dict, Any

def get_tool_detail_str(tool: Dict[str, Any]) -> str:
    """
    A function to get the tool details as a string.
    """
    tool_str = f"Tool Name: {tool['title']}\n"
    if len(tool['link']) > 0:
        tool_str += f"Tool Link: /{tool['link']}\n" if tool['link'][0] != "/" else f"Tool Link: {tool['link']}\n"
    else:
        tool_str += f"Tool Link: /{tool['link']}\n"
    tool_str += f"Tool Description: {tool['description']}\n"
    return tool_str

def format_tool_config_to_prompt(tool_config: List[Dict[str, Any]]) -> str:
    """
    A function to format the tool config to a prompt.
    """
    prompt = """
    Following are the tools that are available for use:

    *Remember to return only tools that are functional not the master tools which have sub-tools. Always be precise about the links as they will be causing redirections*
    
    """ 
    tool_prog = 1
    for t in tool_config:
        prompt += f"{tool_prog}.\n"
        tool_prog += 1
        prompt += get_tool_detail_str(t)
        prompt += "\n"
        if "subElements" in t:
            for s in t["subElements"]:
                prompt += f"{tool_prog}.\n"
                tool_prog += 1
                prompt += get_tool_detail_str(s)
                prompt += "\n"
                if "subElements" in s:
                    for ss in s["subElements"]:
                        prompt += f"{tool_prog}.\n"
                        tool_prog += 1
                        prompt += get_tool_detail_str(ss)
                        prompt += "\n"
    return prompt