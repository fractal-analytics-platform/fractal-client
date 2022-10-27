import asyncio
from copy import deepcopy
from typing import List

from fastapi import APIRouter
from fastapi import Depends
from fastapi import HTTPException
from fastapi import status
from sqlmodel import select

from ...db import AsyncSession
from ...db import DBSyncSession
from ...db import get_db
from ...db import get_sync_db
from ...models import Task
from ...models import TaskCreate
from ...models import TaskRead
from ...models import TaskUpdate
from ...security import current_active_user
from ...security import User

router = APIRouter()


@router.post("/collect/", status_code=status.HTTP_201_CREATED)
async def collect_core_tasks(
    user: User = Depends(current_active_user),
):
    raise NotImplementedError


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


@router.post("/", response_model=TaskRead, status_code=status.HTTP_201_CREATED)
async def create_task(
    task: TaskCreate,
    user: User = Depends(current_active_user),
    db: AsyncSession = Depends(get_db),
):

    # Check that there is no task with the same user and name
    stm = select(Task).where(Task.name == task.name)
    res = await db.execute(stm)
    if res.scalars().all():
        raise HTTPException(
            status_code=status.HTTP_422_UNPROCESSABLE_ENTITY,
            detail=f"Task name ({task.name}) already in use",
        )

    db_task = Task.from_orm(task)
    db.add(db_task)
    await db.commit()
    await db.refresh(db_task)
    return db_task


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
