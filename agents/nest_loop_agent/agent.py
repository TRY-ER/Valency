from google.adk.agents import LoopAgent, Agent
from google.adk.tools import ToolContext
from dotenv import load_dotenv
load_dotenv("../.env")
import os

def increment_tool(x: int) -> int:
    """A simple tool that increments a number by 1."""
    return x + 1

def exit_loop(tool_context: ToolContext) -> dict:
    """A tool that exits the loop."""
    tool_context.actions.escalate = True
    return {}

GEMINI_MODEL = os.getenv("GEMINI_MODEL", "gemini-2.5-flash-preview-05-20")

increment_agent_instruction = """
you always return a number. If you are given an number you increment it by 1 and return the incremented value.
Never return a string or your thoughts and understandings, always return a number.

Always look into the context for the latest incremented value and if current query does not contain any number understand to navigate to the last value and continue the process.

"""

from pydantic import BaseModel

increment_agent = Agent(
    name="IncrementAgent",
    model=GEMINI_MODEL,
    instruction=increment_agent_instruction,
    tools=[increment_tool],
    output_key="incremented_value",
)

nest_review_agent_instruction = """
Review the incremented value and decide whether to exit the loop or continue.
Incremented value = {{incremented_value}}
If the incremented value is greater than 3, exit the loop by calling the exit_loop tool.
Else continue to loop process!
"""

nest_review_agent = Agent(
    name="NestReviewAgent",
    model=GEMINI_MODEL,
    instruction=nest_review_agent_instruction,
    tools=[exit_loop],
)

nest_loop_agent = LoopAgent(
    name="NestedLoopAgent",
    max_iterations=10,
    sub_agents=[increment_agent, nest_review_agent],
    description="A nested loop agent that increments a number until the nest review agent decides to exit.",
)

root_review_agent_instruction = """
Review the findings and decide whether to exit the loop or continue based on the increment value provided.

Incremented value = {{incremented_value}}

If the incremented value is greater than 5 exit the loop by calling the exit_loop tool.
Else continue to loop process !
"""

root_review_agent = Agent(
    name="RootReviewAgent",
    model=GEMINI_MODEL,
    instruction=root_review_agent_instruction,
    tools=[exit_loop],
)


root_agent = LoopAgent(
    name = "RootLoopAgent",
    max_iterations = 100,
    sub_agents = [nest_loop_agent, root_review_agent],
    description = "Root agent that manages the nested loop and review process. The nested agents increase the number given by the user till root review agent decides to exit the loop.",
)