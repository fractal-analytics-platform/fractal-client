from fastapi import APIRouter


router_default = APIRouter()
router_v1 = APIRouter()


@router_default.get("/alive/")
async def alive():
    return dict(alive=True)
