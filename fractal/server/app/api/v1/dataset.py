from fastapi import APIRouter


router = APIRouter()


@router.get("/")
async def get_list_dataset():
    raise NotImplementedError
