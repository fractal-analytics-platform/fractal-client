from fastapi import APIRouter
from ...models.project import ProjectCreateScm


router = APIRouter()


@router.post("/")
async def create_project(
    project: ProjectCreateScm
):
    raise NotImplementedError
