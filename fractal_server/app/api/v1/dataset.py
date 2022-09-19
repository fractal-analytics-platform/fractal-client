from fastapi import APIRouter
from fastapi import Depends
from sqlmodel import select

from ...db import AsyncSession
from ...db import get_db
from ...models import Dataset
from ...models import DatasetRead
from ...models import Project
from ...security import current_active_user
from ...security import User


router = APIRouter()


@router.get("/{project_id}/{dataset_id}", response_model=DatasetRead)
async def get_dataset(
    project_id: int,
    dataset_id: int,
    user: User = Depends(current_active_user),
    db: AsyncSession = Depends(get_db),
):
    stm = (
        select(Project, Dataset)
        .join(Dataset)
        .where(Project.user_owner_id == user.id)
        .where(Project.id == project_id)
        .where(Dataset.id == dataset_id)
    )
    project, dataset = (await db.execute(stm)).one()
    return dataset
