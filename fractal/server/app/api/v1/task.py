import asyncio
from collections import Counter
from typing import Any
from typing import Dict

from fastapi import APIRouter
from fastapi import Depends
from fastapi import status
from sqlalchemy.exc import IntegrityError
from sqlmodel import select

from ....tasks import collect_tasks
from ...db import AsyncSession
from ...db import get_db
from ...models.task import Task
from ...models.task import TaskCreate
from ...security import current_active_user
from ...security import User


router = APIRouter()


async def upsert_task(
    task: Dict[str, Any],
    db: AsyncSession = Depends(get_db),
) -> str:
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


async def collect_tasks_headless(
    db: AsyncSession = Depends(get_db),
):
    out = dict(inserted=0, updated=0)
    results = await asyncio.gather(
        *[upsert_task(task, db) for task in collect_tasks()]
    )
    out.update(dict(Counter(results)))
    return out


@router.post("/collect/", status_code=status.HTTP_201_CREATED)
async def collect_core_tasks(
    user: User = Depends(current_active_user),
    db: AsyncSession = Depends(get_db),
):
    return await collect_tasks_headless(db)


@router.post("/")
async def create_task():
    raise NotImplementedError
