'''
Utility module for reading files.
'''
import os

def read_markdown_file(file_name: str) -> str:
    """
    Reads a markdown file from the 'repo' directory within this package
    and returns its content as a string.

    Args:
        file_name: The name of the markdown file (e.g., "agent_instructions.md")
                   located in the 'utils/repo/' directory.

    Returns:
        The content of the markdown file as a string.
        Returns an error message string if the file cannot be read.
    """
    # Get the directory of the current file (file_reader.py)
    current_dir = os.path.dirname(os.path.abspath(__file__))
    # Construct the path to the 'repo' directory and then to the file
    file_path = os.path.join(current_dir, "repo", file_name)

    try:
        with open(file_path, "r", encoding="utf-8") as f:
            content = f.read()
        return content
    except FileNotFoundError:
        # Try to provide a more helpful path in the error message
        expected_path = os.path.join(os.path.basename(current_dir), "repo", file_name)
        return f"Error: File '{file_name}' not found. Looked for: {expected_path} relative to utils."
    except Exception as e:
        return f"Error reading file '{file_name}': {e}"
