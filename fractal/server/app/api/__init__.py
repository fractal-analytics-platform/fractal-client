from fastapi import APIRouter

from ...config import settings
from .v1.security import router as security_router


router_default = APIRouter()
router_v1 = APIRouter()

router_v1.include_router(security_router, prefix="/auth")


@router_default.get("/alive/")
async def alive():
    return dict(
        alive=True,
        deployment_type=settings.DEPLOYMENT_TYPE,
        version=settings.PROJECT_VERSION,
    )
