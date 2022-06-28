from fastapi import APIRouter

from ...config import settings
from .v1.project import router as project_router


router_default = APIRouter()
router_v1 = APIRouter()

router_v1.include_router(project_router, prefix="/project", tags=["project"])


@router_default.get("/alive/")
async def alive():
    return dict(
        alive=True,
        deployment_type=settings.DEPLOYMENT_TYPE,
        version=settings.PROJECT_VERSION,
    )
