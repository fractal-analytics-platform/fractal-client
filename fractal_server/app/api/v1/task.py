import asyncio
from copy import deepcopy
from pathlib import Path
from typing import List

from fastapi import APIRouter
from fastapi import BackgroundTasks
from fastapi import Depends
from fastapi import status
from sqlmodel import select

from ....config import get_settings
from ....syringe import Inject
from ....tasks.collection import _TaskCollectPip
from ....tasks.collection import create_package_dir_pip
from ....tasks.collection import create_package_environment_pip
from ...db import AsyncSession
from ...db import DBSyncSession
from ...db import get_db
from ...db import get_sync_db
from ...models import Task
from ...schemas import TaskCollectPip
from ...schemas import TaskCreate
from ...schemas import TaskRead
from ...schemas import TaskUpdate
from ...security import current_active_user
from ...security import User

router = APIRouter()


@router.get("/", response_model=List[TaskRead])
async def get_list_task(
    user: User = Depends(current_active_user),
    db: AsyncSession = Depends(get_db),
):
    stm = select(Task)
    res = await db.execute(stm)
    task_list = res.scalars().unique().fetchall()
    await asyncio.gather(*[db.refresh(t) for t in task_list])
    return task_list


@router.get("/{task_id}", response_model=TaskRead)
def get_task(
    task_id: int,
    user: User = Depends(current_active_user),
    db_sync: DBSyncSession = Depends(get_sync_db),
):
    task = db_sync.get(Task, task_id)
    return task


async def _background_collect_pip(
    venv_path: Path, task_pkg: _TaskCollectPip, db: AsyncSession
) -> List[Task]:
    task_list = await create_package_environment_pip(
        venv_path=venv_path, task_pkg=task_pkg
    )
    tasks = await _insert_tasks(task_list=task_list, db=db)
    return tasks


async def _insert_tasks(
    task_list: List[TaskCreate],
    db: AsyncSession,
) -> List[Task]:
    """
    Insert tasks into database
    """
    task_db_list = [Task.from_orm(t) for t in task_list]
    db.add_all(task_db_list)
    await db.commit()
    await asyncio.gather(*[db.refresh(t) for t in task_db_list])
    return task_db_list


@router.post(
    "/pip/",
    status_code=status.HTTP_201_CREATED,
)
async def collect_tasks_pip(
    task_collect: TaskCollectPip,
    background_tasks: BackgroundTasks,
    global_task: bool = True,
    user: User = Depends(current_active_user),
    db: AsyncSession = Depends(get_db),
) -> Task:

    task_pkg = _TaskCollectPip(**task_collect.dict())
    venv_path = create_package_dir_pip(task_pkg=task_pkg)

    background_tasks.add_task(
        _background_collect_pip, venv_path=venv_path, task_pkg=task_pkg, db=db
    )
    settings = Inject(get_settings)
    return dict(
        venv_path=venv_path.relative_to(settings.FRACTAL_ROOT),
        message=f"Installing {task_pkg.package} in the background",
    )


@router.patch("/{task_id}", response_model=TaskRead)
async def patch_task(
    task_id: int,
    task_update: TaskUpdate,
    user: User = Depends(current_active_user),
    db: AsyncSession = Depends(get_db),
):

    # FIXME add user-owned tasks

    db_task = await db.get(Task, task_id)

    for key, value in task_update.dict(exclude_unset=True).items():
        if key == "name":
            setattr(db_task, key, value)
        elif key == "default_args":
            current_default_args = deepcopy(db_task._arguments)
            current_default_args.update(value)
            setattr(db_task, key, current_default_args)
        else:
            raise Exception("patch_task endpoint cannot set {key=}")

    await db.commit()
    await db.refresh(db_task)
    return db_task
