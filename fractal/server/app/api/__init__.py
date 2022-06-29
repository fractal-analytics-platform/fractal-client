from fastapi import APIRouter

from ...config import settings
from .v1.dataset import router as dataset_router
from .v1.project import router as project_router
from .v1.task import router as task_router


router_default = APIRouter()
router_v1 = APIRouter()

router_v1.include_router(project_router, prefix="/project", tags=["project"])
router_v1.include_router(dataset_router, prefix="/dataset", tags=["dataset"])
router_v1.include_router(task_router, prefix="/task", tags=["task"])
router_v1.include_router(task_router, prefix="/workflow", tags=["workflow"])


@router_default.get("/alive/")
async def alive():
    return dict(
        alive=True,
        deployment_type=settings.DEPLOYMENT_TYPE,
        version=settings.PROJECT_VERSION,
    )
