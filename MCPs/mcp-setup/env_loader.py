import os

def load_env_vars(env_path=".env"):
    """Loads environment variables from a .env file line by line."""
    try:
        with open(env_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line and '=' in line and not line.startswith('#'):
                    key_part, value_part_from_split = line.split('=', 1)
                    key = key_part.strip()

                    # value_part_from_split is everything after the first '='.
                    # It needs to be stripped and then parsed for actual value vs. comment.
                    value_str_to_parse = value_part_from_split.strip()

                    # Find the start of a comment (#) not inside quotes
                    comment_idx = -1
                    is_in_quotes = False
                    quote_char_for_iter = ''
                    for i, char_iter in enumerate(value_str_to_parse):
                        if char_iter == "'" or char_iter == '"':
                            if not is_in_quotes:
                                is_in_quotes = True
                                quote_char_for_iter = char_iter
                            elif char_iter == quote_char_for_iter: # Closing quote
                                is_in_quotes = False
                        elif char_iter == '#' and not is_in_quotes: # Comment char outside quotes
                            comment_idx = i
                            break
                    
                    # Extract the value before the comment
                    if comment_idx != -1:
                        # Take the part before the comment and strip whitespace again
                        val_before_comment = value_str_to_parse[:comment_idx].strip()
                    else:
                        # No comment found, the whole (already stripped) string is the value
                        val_before_comment = value_str_to_parse
                    
                    # Remove surrounding quotes from val_before_comment
                    final_env_value = val_before_comment
                    if len(final_env_value) >= 2 and \
                       ((final_env_value.startswith("'") and final_env_value.endswith("'")) or \
                        (final_env_value.startswith('"') and final_env_value.endswith('"'))):
                        final_env_value = final_env_value[1:-1]
                    
                    # Set environment variable only if key is not empty
                    if key:
                        os.environ[key] = final_env_value
    except FileNotFoundError:
        print(f"Warning: {env_path} file not found. Skipping environment variable loading from file.")
    except Exception as e:
        print(f"Error loading {env_path}: {e}")
