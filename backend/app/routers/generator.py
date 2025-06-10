import io
import base64
import random
import os
import httpx
from typing import Dict, List, Generator, AsyncGenerator, Union
from fastapi import APIRouter, File, UploadFile, HTTPException
from fastapi.responses import StreamingResponse, JSONResponse
from engine.generators.BRICSGenerator import BRICSGenerator
import uuid

router = APIRouter()

# Configuration for external ADMET API
ADMET_API_BASE_URL = os.getenv("ADMET_API_BASE_URL", "http://localhost:8053")
ADMET_API_KEY = os.getenv("ADMET_API_KEY", "")

if not ADMET_API_KEY:
    print("WARNING: ADMET_API_KEY not found in environment variables. LSTM generation will not work.")

# In-memory storage for uploaded files and generation parameters
file_storage: Dict[str, Union[List[str], Dict]] = {}

def event_stream(data_lines, gen: BRICSGenerator, is_polymer: bool) -> Generator[str, None, None]:
        for event in gen.stream_generate(data_lines, is_polymer):
            if event["type"] != "completed":
                yield f"data: {event}\n\n"
            else:
                final_data, total_count = event["data"]
                # Write final data to in-memory file
                output_buffer = io.StringIO()
                output_buffer.write("\n".join(final_data))
                output_text = output_buffer.getvalue()
                output_buffer.close()

                # Base64-encode the text file
                encoded_file = base64.b64encode(output_text.encode()).decode()
                returnable = {
                    "file": encoded_file
                }

                # Send total count
                yield f'data: {{\"total\":{total_count}}}\n\n'
                # Send final file
                yield f'data: {returnable}\n\n'

async def generator_stream(content: dict) -> AsyncGenerator[str, None]:
    """
    Stream LSTM generation results from external ADMET API service.
    """
    try:
        # First, set the generation parameters on the external service
        async with httpx.AsyncClient() as client:
            headers = {"Authorization": f"Bearer {ADMET_API_KEY}"}
            
            # Set generation parameters
            set_response = await client.post(
                f"{ADMET_API_BASE_URL}/lstm/set",
                json={
                    "num_gen": content["num_gen"],
                    "input_type": content["input_type"]
                },
                headers=headers
            )
            
            if set_response.status_code != 200:
                error_response = {
                    "type": "error",
                    "data": f"Failed to set generation parameters: {set_response.text}"
                }
                yield f"data: {error_response}\n\n"
                return
            
            set_data = set_response.json()
            gen_id = set_data["id"]
            
            # Stream generation results
            async with client.stream(
                "GET",
                f"{ADMET_API_BASE_URL}/lstm/stream/{gen_id}",
                headers=headers
            ) as stream_response:
                if stream_response.status_code != 200:
                    error_response = {
                        "type": "error", 
                        "data": f"Failed to stream generation: {stream_response.status_code}"
                    }
                    yield f"data: {error_response}\n\n"
                    return
                
                async for chunk in stream_response.aiter_text():
                    if chunk.strip():
                        yield chunk
                        
    except Exception as e:
        error_response = {
            "type": "error",
            "data": str(e)
        }
        yield f"data: {error_response}\n\n"

@router.post("/brics/smiles/upload")
async def upload_smiles(file: UploadFile = File(...)):
    try:
        content = file.file.read().decode("utf-8").splitlines()
        file.file.close()
        
        # Generate unique ID
        file_id = str(uuid.uuid4())
        
        # Store content
        file_storage[file_id] = content
        
        return JSONResponse({
            "status": "success",
            "file_id": file_id,
            "line_count": len(content)
        })
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))

@router.post("/brics/psmiles/upload")
async def upload_psmiles(file: UploadFile = File(...)):
    try:
        content = file.file.read().decode("utf-8").splitlines()
        file.file.close()
        
        # Generate unique ID
        file_id = str(uuid.uuid4())

        print("generated file id:", file_id)
        
        # Store content
        file_storage[file_id] = content
        
        return JSONResponse({
            "status": "success",
            "file_id": file_id,
            "line_count": len(content)
        })
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))


@router.post("/lstm/set")
async def set_generation_psmiles(data: dict):
    try:
        num_gen = data["num_gen"]
        input_type = data["input_type"] 
        gen_id = str(uuid.uuid4())

        file_storage[gen_id] = {
            "input_type": input_type,
            "num_gen": num_gen
        } 

        return JSONResponse({
            "status": "success",
            "id": gen_id,
        })
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))

@router.get("/brics/smiles/stream/{file_id}")
async def stream_smiles(file_id: str):
    if file_id not in file_storage:
        raise HTTPException(status_code=404, detail="File not found")
    
    print("recieved file_id:", file_id)
    content = file_storage[file_id]
    del file_storage[file_id]
    generator = BRICSGenerator()
    return StreamingResponse(
        event_stream(content, generator, is_polymer=False),
        media_type="text/event-stream"
    )

@router.get("/brics/psmiles/stream/{file_id}")
async def stream_psmiles_brics(file_id: str):
    if file_id not in file_storage:
        raise HTTPException(status_code=404, detail="File not found")
    
    content = file_storage[file_id]
    del file_storage[file_id]
    generator = BRICSGenerator()
    return StreamingResponse(
        event_stream(content, generator, is_polymer=True),
        media_type="text/event-stream"
    )

@router.get("/lstm/stream/{gen_id}")
async def stream_psmiles_lstm(gen_id: str):
    if gen_id not in file_storage:
        raise HTTPException(status_code=404, detail="Generation ID not found")
    
    content = file_storage[gen_id]
    del file_storage[gen_id]
    
    # Ensure content is a dict (generation parameters)
    if not isinstance(content, dict):
        raise HTTPException(status_code=400, detail="Invalid generation parameters")
    
    return StreamingResponse(
        generator_stream(content),
        media_type="text/event-stream"
    )

if __name__ == "__main__":
    # Test the external LSTM API connectivity
    import asyncio
    
    async def test_api():
        if not ADMET_API_KEY:
            print("ERROR: ADMET_API_KEY not set in environment variables")
            return
            
        try:
            async with httpx.AsyncClient() as client:
                headers = {"Authorization": f"Bearer {ADMET_API_KEY}"}
                response = await client.get(f"{ADMET_API_BASE_URL}/health", headers=headers)
                if response.status_code == 200:
                    print(f"✓ External ADMET API is healthy: {response.json()}")
                else:
                    print(f"✗ External ADMET API health check failed: {response.status_code}")
        except Exception as e:
            print(f"✗ Failed to connect to external ADMET API: {e}")
    
    asyncio.run(test_api())