from typing import List

from fastapi import APIRouter
from fastapi import Depends
from fastapi import HTTPException
from fastapi import status
from sqlmodel import select

from ...db import AsyncSession
from ...db import get_db
from ...models.project import Dataset
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
    stm = select(Project).where(Project.user_owner_id == user.id)
    res = await db.execute(stm)
    project_list = res.scalars().all()
    return project_list


@router.get("/{project_id}", response_model=ProjectRead)
async def get_project(
    project_id: int,
    user: User = Depends(current_active_user),
    db: AsyncSession = Depends(get_db),
):
    project = await db.get(Project, project_id)
    if not project:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND, detail="Project not found"
        )
    if project.user_owner_id != user.id:
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Not allowed on project",
        )
    return project


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
    db_project.dataset_list.append(Dataset(name=project.default_dataset_name))

    db_project.user_owner_id = user.id
    db.add(db_project)
    await db.commit()
    await db.refresh(db_project)
    return db_project


@router.post("/{project_id}")
async def add_dataset(
    project_id: int,
):
    """
    Add new dataset to current project
    """
    raise NotImplementedError


@router.post("/{project_id}/{dataset_id}")
async def add_resource(
    project_id: str,
    dataset_id: int,
    user: User = Depends(current_active_user),
    db: AsyncSession = Depends(get_db),
):
    raise NotImplementedError


@router.patch("/{project_id}/{dataset_id}")
async def modify_dataset(
    project_id: int,
    dataset_id: int,
):
    raise NotImplementedError


@router.post("/apply/{project_id}/{dataset_id}")
async def apply_workflow():
    raise NotImplementedError
