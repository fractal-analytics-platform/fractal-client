from fastapi import APIRouter
from fastapi import Depends

from ...db import AsyncSession
from ...db import get_db
from ...models.project import Project
from ...security import current_active_user
from ...security import User


router = APIRouter()


@router.post("/", response_model=Project, status_code=201)
async def create_project(
    project: Project,
    user: User = Depends(current_active_user),
    db: AsyncSession = Depends(get_db),
):
    from devtools import debug

    debug(user)
    debug(project)

    db.add(project)
    await db.commit()

    return project
