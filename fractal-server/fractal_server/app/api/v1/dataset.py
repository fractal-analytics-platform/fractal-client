import asyncio

from fastapi import APIRouter
from fastapi import Depends
from fastapi import HTTPException
from fastapi import status
from pydantic import UUID4

from ...db import AsyncSession
from ...db import get_db
from ...models import Dataset
from ...models import DatasetRead
from ...models import LinkUserProject
from ...models import Project
from ...security import current_active_user
from ...security import User


router = APIRouter()


async def get_dataset_check_owner(
    *,
    project_id: int,
    dataset_id: int,
    user_id: UUID4,
    db: AsyncSession = Depends(get_db),
) -> Project:
    """
    Check that user is a member of project and return

    Raise 403_FORBIDDEN if the user is not a member
    Raise 404_NOT_FOUND if the project does not exist
    """
    project, dataset, link_user_project = await asyncio.gather(
        db.get(Project, project_id),
        db.get(Dataset, dataset_id),
        db.get(LinkUserProject, (project_id, user_id)),
    )
    if not project:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND, detail="Project not found"
        )
    if not link_user_project:
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail=f"Not allowed on project {project_id}",
        )
    return dataset


@router.get("/{project_id}/{dataset_id}", response_model=DatasetRead)
async def get_dataset(
    project_id: int,
    dataset_id: int,
    user: User = Depends(current_active_user),
    db: AsyncSession = Depends(get_db),
):

    dataset = await get_dataset_check_owner(
        project_id=project_id, dataset_id=dataset_id, user_id=user.id, db=db
    )

    return dataset
