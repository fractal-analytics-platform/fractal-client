from fastapi import APIRouter

from .v1.security import router as security_router


router_default = APIRouter()
router_v1 = APIRouter()

router_v1.include_router(security_router)


@router_default.get("/alive/")
async def alive():
    return dict(alive=True)
