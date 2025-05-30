from google.adk.tools.mcp_tool.mcp_toolset import (
    MCPToolset,
    SseServerParams
)

async def get_tools_async():
    toolset = MCPToolset(connection_params=SseServerParams(url="http://localhost:8058/sse"))
    tools = await toolset.get_tools()
    print("Tools loaded successfully")
    print("Available tools:", [tool._mcp_tool.name for tool in tools])
    # await toolset.close()
    return tools, toolset.close
    # return tools, exit_stack

async def main():
    tools, closer = await get_tools_async()
    for tool in tools:
        if tool._mcp_tool.name == "get_brics_candidates":
            print("Tool found:", tool._mcp_tool.name)
            response = await tool.run_async(args={"smiles_list": ["CCO", "CCc1ccccc1Br"]}, tool_context=None)
            print("Response from tool:", response)
    await closer()

if __name__ == "__main__":
    import asyncio
    # tools, closer = asyncio.run(get_tools_async())
    asyncio.run(main())