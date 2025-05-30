import os
import sys
import json
from env_loader import load_env_vars

load_env_vars()

from fastmcp import FastMCP
from rxn4chemistry import RXN4ChemistryWrapper

# --- Initialize RXN4Chemistry Wrapper ---
API_KEY = os.getenv("RXN_API_KEY")
BASE_URL = os.getenv("RXN_BASE_URL")
PROJECT_ID = os.getenv("RXN_PROJECT_ID")

if not API_KEY:
    print("Error: RXN_API_KEY environment variable not set.")
    sys.exit(1)

try:
    rxn_wrapper = RXN4ChemistryWrapper(api_key=API_KEY, project_id=PROJECT_ID)
    if BASE_URL:
        rxn_wrapper.set_base_url(BASE_URL)
        print(f"RXN4ChemistryWrapper initialized. Configured to use Base URL: {BASE_URL}")
    else:
        print(f"RXN4ChemistryWrapper initialized with default Base URL.")
    # The following project initialization block was commented out and relied on list_projects,
    # which seems unavailable. It will remain commented.
    # try:
    #     projects = rxn_wrapper.list_projects() # This method causes an error
    #     if projects and projects['content']:
    #         print(f"RXN Wrapper initialized. Available projects found. Set project if needed for specific tools.")
    #     else:
    #         print("RXN Wrapper initialized. No existing projects found or default project created.")
    # except Exception as e:
    #     print(f"Warning: Could not list or set/create default project: {e}")

except Exception as e:
    print(f"Error: RXN API Key is invalid or authentication failed. Please check RXN_API_KEY or other initialization issue: {e}")
    sys.exit(1)

# --- Create MCP Server ---
rxn_mcp_host = os.getenv("RXN_MCP_HOST", "0.0.0.0")
rxn_mcp_port = int(os.getenv("RXN_MCP_PORT", "8053")) # Choose a different port

mcp = FastMCP(
    "RXN for Chemistry MCP Server",
    description="Provides tools to interact with the IBM RXN for Chemistry API.",
    # host=rxn_mcp_host, # Host and port are for run() method
    # port=rxn_mcp_port
)

# --- RXN Tools ---

@mcp.tool()
def predict_reaction_smiles(precursors_smiles: str) -> str:
    """
    Predicts the product of a chemical reaction given the starting materials (precursors) as SMILES.
    Args:
        precursors_smiles: A string containing the SMILES of the precursor molecules, separated by a dot (e.g., 'CCO.Cl').
    Returns:
        A JSON string containing the prediction ID. Use 'get_prediction_results' to get the actual outcome.
    """
    try:
        response = rxn_wrapper.predict_reaction(precursors_smiles)
        return json.dumps({"data": response})
    except Exception as e:
        return json.dumps({"error": f"Failed to predict reaction for '{precursors_smiles}'.", "details": str(e)})

@mcp.tool()
def get_reaction_prediction_results(prediction_id: str) -> str:
    """
    Retrieves the results of a reaction prediction task.
    Args:
        prediction_id: The ID of the prediction task.
    Returns:
        A JSON string containing the prediction results, including status and outcome.
    """
    try:
        results = rxn_wrapper.get_predict_reaction_results(prediction_id)
        return json.dumps({"data": results})
    except Exception as e:
        return json.dumps({"error": f"Failed to get reaction prediction results for ID '{prediction_id}'.", "details": str(e)})

@mcp.tool()
def predict_automatic_retrosynthesis_smiles(product_smiles: str, availability_pricing_threshold: float = 0.0,
                                            available_smiles_list: list[str] | None = None,
                                            exclude_smiles_list: list[str] | None = None,
                                            exclude_substructures_smiles_list: list[str] | None = None,
                                            fap_value: float = 0.6, max_steps: int = 3,
                                            n_beams: int = 3, reaction_class_priority_list: list[str] | None = None,
                                            project_id: str | None = None) -> str:
    """
    Predicts possible retrosynthetic routes for a given target molecule SMILES.
    Args:
        product_smiles: SMILES string of the target molecule.
        availability_pricing_threshold: Exclude precursors with price higher than this (0.0 means no filtering).
        available_smiles_list: List of SMILES that are considered available.
        exclude_smiles_list: List of SMILES to exclude from precursors.
        exclude_substructures_smiles_list: List of SMILES substructures to exclude from precursors.
        fap_value: Forward Acceptance Probability threshold (0.0 to 1.0).
        max_steps: Maximum number of retrosynthetic steps.
        n_beams: Number of beams to use in the search.
        reaction_class_priority_list: List of reaction classes to prioritize.
        project_id: Optional ID of the project to associate this prediction with.
    Returns:
        A JSON string containing the prediction ID. Use 'get_retrosynthesis_prediction_results' to get the routes.
    """
    try:
        # The wrapper might handle None for lists gracefully, but being explicit or converting to empty lists if API expects it.
        # For now, passing them as is.
        current_project_id = project_id or rxn_wrapper.project_id # Use provided or globally set
        if not current_project_id:
             return json.dumps({"error": "Project ID is required for retrosynthesis prediction. Either set a global project or provide 'project_id'.", "details": "No project context."})

        # Temporarily set project for this call if provided and different from current
        original_project_id = rxn_wrapper.project_id
        if project_id and project_id != original_project_id:
            rxn_wrapper.set_project(project_id)

        response = rxn_wrapper.predict_automatic_retrosynthesis(
            product=product_smiles,
            availability_pricing_threshold=availability_pricing_threshold,
            available_smiles_list=available_smiles_list or [], # API might prefer empty list over None
            exclude_smiles_list=exclude_smiles_list or [],
            exclude_substructures_smiles_list=exclude_substructures_smiles_list or [],
            fap_value=fap_value,
            max_steps=max_steps,
            n_beams=n_beams,
            reaction_class_priority_list=reaction_class_priority_list or []
        )
        
        if project_id and project_id != original_project_id: # Reset to original project if it was changed
             if original_project_id:
                rxn_wrapper.set_project(original_project_id)
             # else: # If there was no original project, what to do? For now, leave it as is.
                 # Consider if unsetting a project is possible/needed.

        return json.dumps({"data": response})
    except Exception as e: # Changed to generic Exception
        return json.dumps({"error": f"RXN Error during retrosynthesis: {e}", "details": str(e)})

@mcp.tool()
def get_retrosynthesis_prediction_results(prediction_id: str, project_id: str | None = None) -> str:
    """
    Retrieves the results of an automatic retrosynthesis prediction task.
    Args:
        prediction_id: The ID of the retrosynthesis prediction task.
        project_id: Optional ID of the project this prediction belongs to. If not provided, uses the globally set project.
    Returns:
        A JSON string containing the retrosynthesis prediction results.
    """
    try:
        original_project_id = rxn_wrapper.project_id
        if project_id and project_id != original_project_id:
            rxn_wrapper.set_project(project_id)
        
        results = rxn_wrapper.get_predict_automatic_retrosynthesis_results(prediction_id)
        
        if project_id and project_id != original_project_id: # Reset to original project
            if original_project_id:
                rxn_wrapper.set_project(original_project_id)

        return json.dumps({"data": results})
    except Exception as e: # Changed to generic Exception
        return json.dumps({"error": f"RXN Error retrieving retrosynthesis results: {e}", "details": str(e)})

@mcp.tool()
def predict_reagents_smiles(starting_material_smiles: str, product_smiles: str, project_id: str | None = None) -> str:
    """
    Predicts the reagents needed to convert a given starting material to a given product.
    Args:
        starting_material_smiles: SMILES string of the starting material.
        product_smiles: SMILES string of the product.
        project_id: Optional ID of the project to associate this prediction with.
    Returns:
        A JSON string containing the prediction ID. Use 'get_reagent_prediction_results' to get the actual reagents.
    """
    try:
        original_project_id = rxn_wrapper.project_id
        current_project_id = project_id or original_project_id
        if not current_project_id:
             return json.dumps({"error": "Project ID is required for reagent prediction. Set a global project or provide 'project_id'.", "details": "No project context."})

        if project_id and project_id != original_project_id:
            rxn_wrapper.set_project(project_id)

        response = rxn_wrapper.predict_reagents(starting_material_smiles, product_smiles)
        
        if project_id and project_id != original_project_id:
            if original_project_id:
                rxn_wrapper.set_project(original_project_id)
        
        return json.dumps({"data": response})
    except Exception as e: # Changed to generic Exception
        return json.dumps({"error": f"RXN Error during reagent prediction: {e}", "details": str(e)})

@mcp.tool()
def get_reagent_prediction_results(prediction_id: str, project_id: str | None = None) -> str:
    """
    Retrieves the results of a reagent prediction task.
    Args:
        prediction_id: The ID of the reagent prediction task.
        project_id: Optional ID of the project this prediction belongs to.
    Returns:
        A JSON string containing the reagent prediction results.
    """
    try:
        original_project_id = rxn_wrapper.project_id
        if project_id and project_id != original_project_id:
            rxn_wrapper.set_project(project_id)

        results = rxn_wrapper.get_predict_reagents_results(prediction_id)

        if project_id and project_id != original_project_id:
            if original_project_id:
                rxn_wrapper.set_project(original_project_id)

        return json.dumps({"data": results})
    except Exception as e: # Changed to generic Exception
        return json.dumps({"error": f"RXN Error retrieving reagent prediction results: {e}", "details": str(e)})

@mcp.tool()
def get_atom_mapping_for_reaction_smiles(reaction_smiles_list: list[str]) -> str:
    """
    Performs atom mapping for a list of chemical reactions provided as SMILES strings.
    Args:
        reaction_smiles_list: A list of reaction SMILES strings (e.g., ['CCO.Cl>>CCCl.O']).
    Returns:
        A JSON string containing the atom mapping results.
    """
    try:
        # The predict_reaction_properties tool is used for atom mapping.
        # It expects a list of reaction SMILES.
        response = rxn_wrapper.predict_reaction_properties(reaction_smiles_list)
        return json.dumps({"data": response})
    except Exception as e:
        return json.dumps({"error": "Failed to get atom mapping.", "details": str(e)})

@mcp.tool()
def translate_paragraph_to_actions(paragraph_text: str) -> str:
    """
    Translates a natural language text description of a chemical procedure into machine-readable actions.
    Args:
        paragraph_text: The text paragraph describing the chemical procedure.
    Returns:
        A JSON string containing the extracted actions.
    """
    try:
        response = rxn_wrapper.paragraph_to_actions(paragraph_text)
        return json.dumps({"data": response})
    except Exception as e: # Changed to generic Exception
        return json.dumps({"error": "RXN Error during paragraph to actions translation.", "details": str(e)})

@mcp.tool()
def digitize_reaction_from_file_id(file_id: str, project_id: str | None = None) -> str:
    """
    Starts the digitization process for a reaction scheme from a previously uploaded file.
    Args:
        file_id: The ID of the file uploaded to RXN (e.g., an image of a reaction scheme).
                 Note: Uploading files needs to be handled separately, this tool assumes file_id is known.
        project_id: Optional ID of the project to associate this digitization with.
    Returns:
        A JSON string containing the response from the digitization request (often includes a task ID or status).
    """
    try:
        original_project_id = rxn_wrapper.project_id
        current_project_id = project_id or original_project_id
        if not current_project_id:
             return json.dumps({"error": "Project ID is required for reaction digitization. Set a global project or provide 'project_id'.", "details": "No project context."})

        if project_id and project_id != original_project_id:
            rxn_wrapper.set_project(project_id)

        response = rxn_wrapper.digitize_reaction(file_id) # This might return results directly or a task ID

        if project_id and project_id != original_project_id:
            if original_project_id:
                rxn_wrapper.set_project(original_project_id)
        
        return json.dumps({"data": response})
    except Exception as e: # Changed to generic Exception
        return json.dumps({"error": f"RXN Error during reaction digitization for file ID '{file_id}'.", "details": str(e)})

# --- File Upload Tool ---

@mcp.tool()
def rxn_upload_file(server_file_path: str) -> str:
    """Uploads a file to the RXN server.
    Args:
        server_file_path: The absolute path to the file on the server where this MCP is running.
                          The file must be accessible by this server process.
    Returns:
        A JSON string with the file ID upon successful upload.
    """
    try:
        # Ensure the file exists before attempting to upload
        if not os.path.exists(server_file_path):
            return json.dumps({"error": "File not found on server.", "details": f"Path '{server_file_path}' does not exist."})
        if not os.path.isfile(server_file_path):
            return json.dumps({"error": "Path is not a file.", "details": f"Path '{server_file_path}' is not a file."})

        response = rxn_wrapper.upload_file(server_file_path)
        return json.dumps({"data": response}) # response usually contains {'response': {'payload': {'id': 'file_id'}}}
    except FileNotFoundError:
        return json.dumps({"error": "File not found during RXN upload process.", "details": f"Path '{server_file_path}' was initially found but disappeared."})
    except Exception as e:
        return json.dumps({"error": f"Failed to upload file '{server_file_path}' to RXN.", "details": str(e)})

# --- Batch Reaction Prediction Tools ---

@mcp.tool()
def predict_reaction_batch_smiles(precursors_smiles_list: list[str]) -> str:
    """
    Predicts products for a batch of chemical reactions given lists of starting materials (precursors) as SMILES.
    Note: Batch predictions are not stored in RXN projects.
    Args:
        precursors_smiles_list: A list of strings, where each string contains SMILES of precursor molecules for one reaction, separated by a dot (e.g., ['CCO.Cl', 'c1ccccc1.BrBr']).
    Returns:
        A JSON string containing the task ID for the batch prediction. Use 'get_reaction_batch_prediction_results' to get outcomes.
    """
    try:
        response = rxn_wrapper.predict_reaction_batch(precursors_list=precursors_smiles_list)
        return json.dumps({"data": response}) # Expected: {'task_id': 'some_id'}
    except Exception as e:
        return json.dumps({"error": "Failed to submit batch reaction prediction.", "details": str(e)})

@mcp.tool()
def get_reaction_batch_prediction_results(task_id: str) -> str:
    """
    Retrieves the results of a batch reaction prediction task.
    Args:
        task_id: The ID of the batch prediction task.
    Returns:
        A JSON string containing the prediction results for the batch.
    """
    try:
        results = rxn_wrapper.get_predict_reaction_batch_results(task_id)
        return json.dumps({"data": results})
    except Exception as e:
        return json.dumps({"error": f"Failed to get batch reaction prediction results for task ID '{task_id}'.", "details": str(e)})

@mcp.tool()
def predict_reaction_batch_topn_smiles(precursors_smiles_lists: list[list[str]], top_n: int) -> str:
    """
    Predicts the top N products for a batch of chemical reactions.
    Note: Batch predictions are not stored in RXN projects.
    Args:
        precursors_smiles_lists: A list of lists of strings. Each inner list contains SMILES of precursor molecules for one reaction (e.g., [['BrBr', 'c1ccccc1'], ['Cl', 'CCO']]).
        top_n: The number of top predictions to return for each reaction.
    Returns:
        A JSON string containing the task ID. Use 'get_reaction_batch_topn_prediction_results' to get outcomes.
    """
    try:
        response = rxn_wrapper.predict_reaction_batch_topn(precursors_lists=precursors_smiles_lists, topn=top_n)
        return json.dumps({"data": response}) # Expected: {'task_id': 'some_id'}
    except Exception as e:
        return json.dumps({"error": "Failed to submit batch top-N reaction prediction.", "details": str(e)})

@mcp.tool()
def get_reaction_batch_topn_prediction_results(task_id: str) -> str:
    """
    Retrieves the results of a batch top-N reaction prediction task.
    Args:
        task_id: The ID of the batch top-N prediction task.
    Returns:
        A JSON string containing the top-N prediction results for the batch.
    """
    try:
        results = rxn_wrapper.get_predict_reaction_batch_topn_results(task_id)
        return json.dumps({"data": results})
    except Exception as e:
        return json.dumps({"error": f"Failed to get batch top-N reaction prediction results for task ID '{task_id}'.", "details": str(e)})

# --- Synthesis Planning Tools ---

@mcp.tool()
def create_synthesis_from_sequence(sequence_id: str, project_id: str | None = None) -> str:
    """
    Creates a synthesis plan in RXN from a retrosynthesis sequence ID.
    Args:
        sequence_id: The ID of the retrosynthetic sequence (obtained from retrosynthesis results).
        project_id: Optional ID of the project to associate this synthesis with. If None, uses the globally set project.
    Returns:
        A JSON string with the details of the created synthesis, including its ID.
    """
    try:
        original_project_id = rxn_wrapper.project_id
        current_project_id = project_id or original_project_id
        if not current_project_id:
            return json.dumps({"error": "Project ID is required for creating a synthesis. Set a global project or provide 'project_id'."})

        if project_id and project_id != original_project_id:
            rxn_wrapper.set_project(project_id)

        response = rxn_wrapper.create_synthesis_from_sequence(sequence_id=sequence_id)
        
        if project_id and project_id != original_project_id:
            if original_project_id: rxn_wrapper.set_project(original_project_id)
        
        return json.dumps({"data": response}) # Expected: {'synthesis_id': 'some_id'}
    except Exception as e: # Changed to generic Exception
        return json.dumps({"error": f"RXN Error during synthesis creation from sequence ID '{sequence_id}'.", "details": str(e)})

@mcp.tool()
def get_synthesis_node_ids(synthesis_id: str, project_id: str | None = None) -> str:
    """
    Retrieves all node IDs for a given synthesis plan.
    Args:
        synthesis_id: The ID of the synthesis plan.
        project_id: Optional ID of the project this synthesis belongs to. If None, uses the globally set project.
    Returns:
        A JSON string containing a list of node IDs.
    """
    try:
        original_project_id = rxn_wrapper.project_id
        current_project_id = project_id or original_project_id
        if not current_project_id:
            return json.dumps({"error": "Project ID is required to get synthesis node IDs. Set a global project or provide 'project_id'."})

        if project_id and project_id != original_project_id:
            rxn_wrapper.set_project(project_id)

        node_ids = rxn_wrapper.get_node_ids(synthesis_id=synthesis_id)
        
        if project_id and project_id != original_project_id:
            if original_project_id: rxn_wrapper.set_project(original_project_id)

        return json.dumps({"data": node_ids})
    except Exception as e: # Changed to generic Exception
        return json.dumps({"error": f"RXN Error retrieving synthesis node IDs for synthesis ID '{synthesis_id}'.", "details": str(e)})

@mcp.tool()
def get_synthesis_node_reaction_settings(synthesis_id: str, node_id: str, project_id: str | None = None) -> str:
    """
    Retrieves the reaction settings (actions and product) for a specific node in a synthesis plan.
    Args:
        synthesis_id: The ID of the synthesis plan.
        node_id: The ID of the node within the synthesis plan.
        project_id: Optional ID of the project. If None, uses the globally set project.
    Returns:
        A JSON string containing the actions and product SMILES for the specified node.
    """
    try:
        original_project_id = rxn_wrapper.project_id
        current_project_id = project_id or original_project_id
        if not current_project_id:
            return json.dumps({"error": "Project ID is required to get node settings. Set a global project or provide 'project_id'."})

        if project_id and project_id != original_project_id:
            rxn_wrapper.set_project(project_id)

        settings = rxn_wrapper.get_reaction_settings(synthesis_id=synthesis_id, node_id=node_id)
        
        if project_id and project_id != original_project_id:
            if original_project_id: rxn_wrapper.set_project(original_project_id)

        return json.dumps({"data": settings}) # Expected: {'actions': [...], 'product': 'smiles'}
    except Exception as e: # Changed to generic Exception
        return json.dumps({"error": f"RXN Error retrieving node settings for node '{node_id}' in synthesis '{synthesis_id}'.", "details": str(e)})

@mcp.tool()
def update_synthesis_node_reaction_settings(synthesis_id: str, node_id: str, actions: list[dict], product_smiles: str, project_id: str | None = None) -> str:
    """
    Updates the reaction settings (actions and product) for a specific node in a synthesis plan.
    Args:
        synthesis_id: The ID of the synthesis plan.
        node_id: The ID of the node within the synthesis plan.
        actions: A list of action dictionaries defining the procedure. Refer to RXN documentation for the structure.
        product_smiles: The SMILES string of the product for this reaction step.
        project_id: Optional ID of the project. If None, uses the globally set project.
    Returns:
        A JSON string confirming the update or an error message.
    """
    try:
        original_project_id = rxn_wrapper.project_id
        current_project_id = project_id or original_project_id
        if not current_project_id:
            return json.dumps({"error": "Project ID is required to update node settings. Set a global project or provide 'project_id'."})

        if project_id and project_id != original_project_id:
            rxn_wrapper.set_project(project_id)

        # The wrapper function is update_reaction_settings
        response = rxn_wrapper.update_reaction_settings(synthesis_id=synthesis_id, node_id=node_id, actions=actions, product=product_smiles)
        
        if project_id and project_id != original_project_id:
            if original_project_id: rxn_wrapper.set_project(original_project_id)

        return json.dumps({"data": response}) # Response structure might vary, often a success message or updated entity
    except Exception as e: # Changed to generic Exception
        return json.dumps({"error": f"RXN Error updating node settings for node '{node_id}' in synthesis '{synthesis_id}'. Check actions format. Details: {str(e)}"})

# --- MCP Server Project Management (Optional, but good practice) ---
@mcp.tool()
def rxn_create_project(name: str) -> str:
    """Creates a new project in RXN for Chemistry.
    Args:
        name: The name for the new project.
    Returns:
        A JSON string with the details of the created project, including its ID.
    """
    try:
        project_details = rxn_wrapper.create_project(name)
        return json.dumps({"data": project_details})
    except Exception as e:
        return json.dumps({"error": f"Failed to create RXN project '{name}'.", "details": str(e)})

@mcp.tool()
def rxn_set_current_project(project_id: str) -> str:
    """Sets the current active project for the RXN wrapper instance for subsequent calls.
    Args:
        project_id: The ID of the project to set as active.
    Returns:
        A JSON string confirming the action or an error.
    """
    try:
        rxn_wrapper.set_project(project_id)
        return json.dumps({"data": f"RXN project successfully set to '{project_id}'. Current project ID: {rxn_wrapper.project_id}"})
    except Exception as e: # Changed to generic Exception
        return json.dumps({"error": f"Failed to set RXN project to '{project_id}'. Details: {str(e)}"})

@mcp.tool()
def rxn_get_current_project_id() -> str:
    """Gets the ID of the currently active project in the RXN wrapper instance.
    Returns:
        A JSON string with the current project ID or a message if no project is set.
    """
    try:
        project_id = rxn_wrapper.project_id
        if project_id:
            return json.dumps({"data": {"current_project_id": project_id}})
        else:
            return json.dumps({"data": {"message": "No project is currently set in the RXN wrapper."}})
    except Exception as e:
        return json.dumps({"error": "Failed to get current RXN project ID.", "details": str(e)})


# --- Main execution ---
if __name__ == "__main__":
    if not API_KEY:
        print("RXN_API_KEY is not set. Server cannot start.", file=sys.stderr)
        sys.exit(1)
    
    print(f"Starting RXN for Chemistry MCP Server...")
    print(f"Server Name: {mcp.name}")
    print(f"Attempting to run on host {rxn_mcp_host}, port {rxn_mcp_port} using SSE transport.")
    sys.stdout.flush()
    try:
        mcp.run(transport="sse", host=rxn_mcp_host, port=rxn_mcp_port)
    except Exception as e:
        print(f"Failed to start MCP server: {e}", file=sys.stderr)
        sys.exit(1)

