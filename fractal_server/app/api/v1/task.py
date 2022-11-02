import asyncio
from copy import deepcopy
from typing import List

from fastapi import APIRouter
from fastapi import Depends
from fastapi import status
from sqlmodel import select

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


async def collect_tasks_headless(
    task_list: List[TaskCreate],
    db: AsyncSession,
    global_task: bool = True,
) -> List[Task]:
    task_db_list = [Task.from_orm(t) for t in task_list]
    db.add_all(task_db_list)
    await db.commit()
    await asyncio.gather(*[db.refresh(t) for t in task_db_list])
    return task_db_list


@router.post(
    "/pip/",
    response_model=List[TaskRead],
    status_code=status.HTTP_201_CREATED,
)
async def collect_tasks_pip(
    task_collect: TaskCollectPip,
    collection_type: str = "pip",
    global_task: bool = True,
    user: User = Depends(current_active_user),
    db: AsyncSession = Depends(get_db),
) -> Task:

    task_pkg = _TaskCollectPip(**task_collect.dict())
    venv_path = create_package_dir_pip(task_pkg=task_pkg)
    task_list = await create_package_environment_pip(
        venv_path=venv_path, task_pkg=task_pkg
    )
    task_db_list = await collect_tasks_headless(
        task_list=task_list,
        db=db,
        global_task=global_task,
    )
    return task_db_list


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
