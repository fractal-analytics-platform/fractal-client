from typing import List

from fastapi import APIRouter
from fastapi import Depends
from sqlalchemy.future import select

from ...db import AsyncSession
from ...db import get_db
from ...models.project import Project
from ...security import current_active_user
from ...security import User


router = APIRouter()


@router.get("/", response_model=List[Project])
async def get_list_project(
    db: AsyncSession = Depends(get_db),
):
    print("I'm here")
    result = await db.execute(select(Project))
    return result.scalars().all()


@router.post("/", response_model=Project, status_code=201)
async def create_project(
    project: Project,
    user: User = Depends(current_active_user),
    db: AsyncSession = Depends(get_db),
):
    project = Project(name="project", project_dir="/tpm/")
    db.add(project)
    await db.commit()
    await db.refresh(project)
    return project
