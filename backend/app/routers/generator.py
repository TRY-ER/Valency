import io
import base64
import random
from typing import Dict, List, Generator
from fastapi import APIRouter, File, UploadFile, HTTPException
from fastapi.responses import StreamingResponse, JSONResponse
from engine.generators.BRICSGenerator import BRICSGenerator
from engine.generators.LSTMGenerator import RNNPolymerGenerator 
import uuid
import torch

router = APIRouter()

# In-memory storage for uploaded files
file_storage: Dict[str, List[str]] = {}

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

def generator_stream(content: dict):
    try:
        if content["input_type"] == "psmiles":
            generator = RNNPolymerGenerator(input_type="psmiles")
            generator.load_model_from_ckpt("./resources/PSMILES_LSTM_1M_5_epochs.pth")
        elif content["input_type"] == "wdg":
            generator = RNNPolymerGenerator(input_type="wdg")
            generator.load_model_from_ckpt("./resources/WDGraph_LSTM_42K_50_epochs.pth")
        for entity in generator.stream_generate(int(content["num_gen"])):
            yield f"data: {entity}\n\n"
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
    
    content = file_storage[file_id]
    del file_storage[file_id]
    generator = BRICSGenerator()
    return StreamingResponse(
        event_stream(content, generator, is_polymer=False),
        media_type="text/event-stream"
    )

@router.get("/brics/psmiles/stream/{file_id}")
async def stream_psmiles(file_id: str):
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
async def stream_psmiles(gen_id: str):
    if gen_id not in file_storage:
        raise HTTPException(status_code=404, detail="File not found")
    
    content = file_storage[gen_id]
    del file_storage[gen_id]
    return StreamingResponse(
        generator_stream(content),
        media_type="text/event-stream"
    )

if __name__ == "__main__":
    generator = RNNPolymerGenerator(input_type="psmiles")
    generator.load_model_from_ckpt("./resources/PSMILES_LSTM_1M_5_epochs.pth")
    for gen in generator.stream_generate(10):
        print('generations >>', gen)