import asyncio
from collections import Counter
from copy import deepcopy
from typing import Any
from typing import Dict
from typing import List

from fastapi import APIRouter
from fastapi import Depends
from fastapi import HTTPException
from fastapi import status
from sqlalchemy.exc import IntegrityError
from sqlmodel import select

from ....tasks import collect_tasks
from ...db import async_session_maker
from ...db import AsyncSession
from ...db import DBSyncSession
from ...db import get_db
from ...db import get_sync_db
from ...models import SubtaskCreate
from ...models import Task
from ...models import TaskCreate
from ...models import TaskRead
from ...models import TaskUpdate
from ...security import current_active_user
from ...security import User

router = APIRouter()


async def upsert_task(
    task: Dict[str, Any],
) -> str:
    async with async_session_maker() as db:
        task_obj = TaskCreate(**task)
        try:
            task_orm = Task.from_orm(task_obj)
            db.add(task_orm)
            await db.commit()
            return "inserted"
        except IntegrityError:
            await db.rollback()
            stm = select(Task).where(Task.name == task_obj.name)
            res = await db.execute(stm)
            this_task = res.scalars().one()
            for key, value in task_obj.dict(exclude={"subtask_list"}).items():
                setattr(this_task, key, value)
            db.add(this_task)
            await db.commit()
            return "updated"


async def collect_tasks_headless():
    out = dict(inserted=0, updated=0)
    results = await asyncio.gather(
        *[upsert_task(task) for task in collect_tasks()]
    )
    results = sorted(results)
    out.update(dict(Counter(results)))
    return out


@router.post("/collect/", status_code=status.HTTP_201_CREATED)
async def collect_core_tasks(
    user: User = Depends(current_active_user),
):
    return await collect_tasks_headless()


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


@router.post(
    "/{parent_task_id}/subtask/",
    response_model=TaskRead,
    status_code=status.HTTP_201_CREATED,
)
async def add_subtask(
    parent_task_id: int,
    subtask: SubtaskCreate,
    user: User = Depends(current_active_user),
    db: AsyncSession = Depends(get_db),
    db_sync: DBSyncSession = Depends(get_sync_db),
):
    parent_task = await db.get(Task, parent_task_id)
    await parent_task.add_subtask(
        db=db,
        **subtask.dict(exclude={"parent_task_id"}),
    )

    # Fetch parent task synchronously to lazily resolve children
    parent_task = db_sync.get(Task, parent_task_id)
    return parent_task
