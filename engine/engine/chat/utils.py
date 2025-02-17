def postproces_response_for_stream(data: str):
    """
    Postprocess the response for stream
    """
    # return data.strip().replace("\n", " ").replace("\r", " ") + "\n\n"
    return data.strip()