import os
import httpx
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel

router = APIRouter()

# Configuration for external ADMET API
ADMET_API_BASE_URL = os.getenv("ADMET_API_BASE_URL", "http://localhost:8053")
ADMET_API_KEY = os.getenv("ADMET_API_KEY", "")

if not ADMET_API_KEY:
    print("WARNING: ADMET_API_KEY not found in environment variables. Similarity search will not work.")

# Define allowed types locally since we removed the engine import
ALLOWED_TYPES = ["MOL", "POLY", "PROT"]

class SimilarityQuery(BaseModel):
    input_type: str
    data: str
    k: int

async def call_remote_similarity_search(payload: SimilarityQuery):
    """
    Call the remote ADMET API similarity search endpoint.
    """
    try:
        async with httpx.AsyncClient() as client:
            headers = {"Authorization": f"Bearer {ADMET_API_KEY}"}

            print("admet api base url >>", ADMET_API_BASE_URL) 
            print("admet key >>", ADMET_API_KEY) 
            response = await client.post(
                f"{ADMET_API_BASE_URL}/ssearch/local",
                json={
                    "input_type": payload.input_type,
                    "data": payload.data,
                    "k": payload.k
                },
                headers=headers
            )
            
            if response.status_code != 200:
                error_detail = f"Remote API error: {response.status_code}"
                try:
                    error_data = response.json()
                    if "detail" in error_data:
                        error_detail = error_data["detail"]
                except:
                    error_detail = f"Remote API error: {response.text}"
                print("error details >>", error_detail)
                raise HTTPException(status_code=response.status_code, detail=error_detail)
            
            return response.json()
            
    except httpx.RequestError as e:
        raise HTTPException(
            status_code=503, 
            detail=f"Failed to connect to remote similarity search service: {str(e)}"
        )
    except httpx.HTTPStatusError as e:
        raise HTTPException(
            status_code=e.response.status_code,
            detail=f"Remote API returned error: {e.response.text}"
        )


@router.post("/ssearch/local")
async def similarity_search(payload: SimilarityQuery):
    """
    Perform similarity search using remote ADMET API service.
    """
    if payload.k < 1:
        raise HTTPException(status_code=400, detail="k should be a positive integer")
    if payload.input_type not in ALLOWED_TYPES:
        raise HTTPException(status_code=400, detail=f"Validation type should be one of {ALLOWED_TYPES}")
    
    # Call the remote similarity search service
    return await call_remote_similarity_search(payload)

if __name__ == "__main__":
    # Test the external ADMET API connectivity
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
    