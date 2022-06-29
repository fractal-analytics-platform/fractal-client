from fastapi import APIRouter


router = APIRouter()


@router.post("/")
async def create_task():
    raise NotImplementedError
