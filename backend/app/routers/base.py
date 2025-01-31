from fastapi import APIRouter

router = APIRouter()

@router.get("/ping")
def get_items():
    return {"response": "working"}