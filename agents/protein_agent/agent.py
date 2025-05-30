# filepath: /home/kalki/src/valency/agents/protein_agent/agent.py
from adk.agent import Agent
from adk.message import Message, MessageType # Added Message import
from sub_agents import (
    alphafold_agent,
    rcsb_agent
)
import asyncio # Added asyncio import

# Define the ProteinAgent
class ProteinAgent(Agent):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        # It's good practice to ensure sub_agents are initialized correctly.
        # Assuming they have default constructors or you'll pass necessary kwargs.
        self.alphafold_agent = alphafold_agent.AlphafoldAgent()
        self.rcsb_agent = rcsb_agent.RCSBAgent()
        # You can load the instructions from the markdown file you created.
        # For example, using a utility function if you have one:
        # self.instructions = load_instructions_from_file("protein_agent_instructions.md")
        # For now, let's print a confirmation or store a simple instruction string.
        self.instructions = "ProteinAgent: I can route queries to AlphaFold or RCSB PDB sub-agents."
        print(f"ProteinAgent initialized with instructions: {self.instructions}")

    async def on_message(self, message: Message, history: list[Message]): # Added type hints
        if message.type == MessageType.USER:
            user_query = message.text.lower()

            # Keywords for AlphaFold/UniProt
            alphafold_keywords = ["alphafold", "uniprot", "prediction", "sequence function", "p0dp23", "q5vsl9", "p12345", "p00766"]
            # Keywords for RCSB/PDB
            rcsb_keywords = ["pdb", "rcsb", "experimental structure", "x-ray", "nmr", "cryo-em", "1a2b", "4hhb", "1tup"]

            use_alphafold = any(keyword in user_query for keyword in alphafold_keywords)
            use_rcsb = any(keyword in user_query for keyword in rcsb_keywords)

            # If specific identifiers are present, prioritize them
            # (This is a simplified check; more robust parsing might be needed)
            if any(pdb_id.lower() in user_query for pdb_id in ["1a2b", "4hhb", "1tup"]) and not use_alphafold:
                use_rcsb = True
                use_alphafold = False
            elif any(uniprot_id.lower() in user_query for uniprot_id in ["p0dp23", "q5vsl9", "p12345", "p00766"]) and not use_rcsb:
                use_alphafold = True
                use_rcsb = False

            response_text = ""

            if use_alphafold and not use_rcsb:
                # Forward to alphafold_agent
                # Assuming sub-agents also have an on_message or similar async method
                # and return a Message object or a string that can be wrapped in a Message.
                alphafold_response_message = await self.alphafold_agent.on_message(message, history) # Pass the original message and history
                if isinstance(alphafold_response_message, Message):
                    response_text = alphafold_response_message.text
                elif isinstance(alphafold_response_message, str):
                    response_text = alphafold_response_message
                else:
                    response_text = "AlphaFold agent did not return a recognizable response."
                await self.send_message(response_text, message_type=MessageType.ASSISTANT)

            elif use_rcsb and not use_alphafold:
                # Forward to rcsb_agent
                rcsb_response_message = await self.rcsb_agent.on_message(message, history) # Pass the original message and history
                if isinstance(rcsb_response_message, Message):
                    response_text = rcsb_response_message.text
                elif isinstance(rcsb_response_message, str):
                    response_text = rcsb_response_message
                else:
                    response_text = "RCSB agent did not return a recognizable response."
                await self.send_message(response_text, message_type=MessageType.ASSISTANT)

            elif use_alphafold and use_rcsb:
                # More complex scenario: potentially needs orchestration
                await self.send_message(
                    "Your query seems to involve both UniProt/AlphaFold and PDB data. Could you please specify which aspect you'd like to focus on first, or I can try to address them sequentially.",
                    message_type=MessageType.ASSISTANT
                )
                # Example: Try AlphaFold first, then PDB based on results (needs more logic)
                # alphafold_response = await self.alphafold_agent.on_message(message, history)
                # await self.send_message(f"From AlphaFold/UniProt: {alphafold_response.text}", message_type=MessageType.ASSISTANT)
                # Potentially use info from alphafold_response to query rcsb_agent
                # rcsb_follow_up_query_text = f"Based on AlphaFold data ({alphafold_response.text}), now query PDB." # Placeholder
                # rcsb_follow_up_message = Message(text=rcsb_follow_up_query_text, type=MessageType.USER) # Create a new message
                # rcsb_response = await self.rcsb_agent.on_message(rcsb_follow_up_message, history)
                # await self.send_message(f"From PDB: {rcsb_response.text}", message_type=MessageType.ASSISTANT)

            else:
                # Default or ambiguous query
                await self.send_message(
                    "I can help with protein information from AlphaFold/UniProt (for predicted structures and sequences) or RCSB PDB (for experimental structures). Please specify what you are looking for, or provide a UniProt ID or PDB ID.",
                    message_type=MessageType.ASSISTANT
                )
        else:
            # Handle other message types if necessary
            pass

if __name__ == '__main__':
    # This is a basic way to run the agent.
    # In a real application, you'd integrate this into your ADK setup.
    agent = ProteinAgent(agent_id="protein_agent_test") # agent_id is often required by base Agent
    print("ProteinAgent initialized. It can use alphafold_agent and rcsb_agent.")
    print("To use it, you would typically send messages to its on_message handler.")

    # Example of how you might simulate a message (requires an event loop)
    async def main():
        # Ensure sub-agents are also initialized if they have async setup needs
        # await agent.alphafold_agent.initialize() # if they have such methods
        # await agent.rcsb_agent.initialize()      # if they have such methods

        sample_message = Message(text="Tell me about UniProt P0DP23", type=MessageType.USER, sender_id="user1", agent_id="protein_agent_test")
        history = []
        await agent.on_message(sample_message, history)

    # Python 3.7+
    # asyncio.run(main())
    # For older versions or different event loop management, adjust accordingly.
    # If running in an environment that already has an asyncio loop (e.g., Jupyter),
    # you might just call `await main()` if in an async context.