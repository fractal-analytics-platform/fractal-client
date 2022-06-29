from typing import List

from fastapi import APIRouter
from fastapi import Depends
from sqlalchemy.future import select

from ...db import AsyncSession
from ...db import get_db
from ...models.project import Project
from ...models.project import ProjectCreate
from ...models.project import ProjectRead
from ...security import current_active_user
from ...security import User


router = APIRouter()


@router.get("/", response_model=List[ProjectRead])
async def get_list_project(
    user: User = Depends(current_active_user),
    db: AsyncSession = Depends(get_db),
):
    result = await db.execute(select(Project))
    return result.scalars().all()


@router.post("/", response_model=ProjectRead, status_code=201)
async def create_project(
    project: ProjectCreate,
    user: User = Depends(current_active_user),
    db: AsyncSession = Depends(get_db),
):
    """
    Create new poject
    """
    db_project = Project.from_orm(project)

    db_project.user_owner_id = user.id
    db.add(db_project)
    await db.commit()
    await db.refresh(db_project)
    return db_project


@router.post("/apply/{project_slug}/{dataset_id}")
async def apply_workflow():
    raise NotImplementedError


@router.post("/{project_slug}")
async def add_dataset():
    """
    Add new dataset to current project
    """
    raise NotImplementedError


@router.post("/{project_slug}/{dataset_id}/")
async def add_resource():
    raise NotImplementedError


@router.patch("/{project_slug}/{dataset_id}")
async def modify_dataset():
    raise NotImplementedError
