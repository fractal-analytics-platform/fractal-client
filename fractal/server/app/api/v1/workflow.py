from fastapi import APIRouter
from fastapi import Depends

from ...db import AsyncSession
from ...db import get_db
from ...security import current_active_user
from ...security import User

router = APIRouter()


@router.post("/")
async def create_workflow(
    user: User = Depends(current_active_user),
    db: AsyncSession = Depends(get_db),
):
    """
    Add new global workflow
    """
    raise NotImplementedError


@router.post("/{wofklow_id}")
async def add_task(
    user: User = Depends(current_active_user),
    db: AsyncSession = Depends(get_db),
):
    """
    Add global task to global workflow
    """
    raise NotImplementedError


@router.get("/")
async def get_list_workflow(
    user: User = Depends(current_active_user),
    db: AsyncSession = Depends(get_db),
):
    """
    List public workflows
    """
    raise NotImplementedError
