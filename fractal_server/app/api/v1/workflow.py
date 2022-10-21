"""
Copyright 2022 (C) Friedrich Miescher Institute for Biomedical Research and
University of Zurich

Original author(s):
Jacopo Nespolo <jacopo.nespolo@exact-lab.it>

This file is part of Fractal and was originally developed by eXact lab S.r.l.
<exact-lab.it> under contract with Liberali Lab from the Friedrich Miescher
Institute for Biomedical Research and Pelkmans Lab from the University of
Zurich.
"""
from fastapi import APIRouter
from fastapi import Depends
from fastapi import HTTPException
from fastapi import Response
from fastapi import status
from sqlmodel import select

from ...db import AsyncSession
from ...db import get_db
from ...models import Workflow
from ...models import WorkflowCreate
from ...models import WorkflowRead
from ...models import WorkflowTaskCreate
from ...models import WorkflowUpdate
from ...security import current_active_user
from ...security import User
from .project import get_project_check_owner

router = APIRouter()


@router.post(
    "/", response_model=WorkflowRead, status_code=status.HTTP_201_CREATED
)
async def create_workflow(
    workflow: WorkflowCreate,
    user: User = Depends(current_active_user),
    db: AsyncSession = Depends(get_db),
):
    await get_project_check_owner(
        project_id=workflow.project_id,
        user_id=user.id,
        db=db,
    )
    # Check that there is no workflow with the same name
    # and same project_id
    stm = (
        select(Workflow)
        .where(Workflow.name == workflow.name)
        .where(Workflow.project_id == workflow.project_id)
    )
    res = await db.execute(stm)
    if res.scalars().all():
        raise HTTPException(
            status_code=status.HTTP_422_UNPROCESSABLE_ENTITY,
            detail=f"Workflow with name={workflow.name} and\
                    project_id={workflow.project_id} already in use",
        )
    db_workflow = Workflow.from_orm(workflow)
    db.add(db_workflow)
    await db.commit()
    await db.refresh(db_workflow)
    return db_workflow


@router.delete("/{_id}/", status_code=status.HTTP_204_NO_CONTENT)
async def delete_workflow(
    _id: int,
    user: User = Depends(current_active_user),
    db: AsyncSession = Depends(get_db),
):
    workflow = await db.get(Workflow, _id)
    if not workflow:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND, detail="Workflow not found"
        )
    await get_project_check_owner(
        project_id=workflow.project_id,
        user_id=user.id,
        db=db,
    )
    await db.delete(workflow)
    await db.commit()

    return Response(status_code=status.HTTP_204_NO_CONTENT)


@router.get("/{_id}/", response_model=WorkflowRead)
async def get_workflow(
    _id: int,
    user: User = Depends(current_active_user),
    db: AsyncSession = Depends(get_db),
):
    # TODO move check autorization as first thing (issue #171)
    workflow = await db.get(Workflow, _id)
    if not workflow:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND, detail="Workflow not found"
        )
    # check autorization
    await get_project_check_owner(
        project_id=workflow.project_id,
        user_id=user.id,
        db=db,
    )
    return workflow


@router.post(
    "/{_id}/add-task/",
    response_model=WorkflowRead,
    status_code=status.HTTP_201_CREATED,
)
async def add_task_to_workflow(
    _id: int,
    new_task: WorkflowTaskCreate,
    user: User = Depends(current_active_user),
    db: AsyncSession = Depends(get_db),
):
    # TODO move check autorization as first thing (issue #171)
    workflow = await db.get(Workflow, _id)
    if not workflow:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND, detail="Workflow not found"
        )
    # check autorization
    await get_project_check_owner(
        project_id=workflow.project_id,
        user_id=user.id,
        db=db,
    )
    await workflow.insert_task(new_task.task_id, db=db)

    db.merge(workflow)
    db.commit()
    await db.refresh(workflow)

    return workflow


@router.delete(
    "{_id}/rm-task/{task_id}", status_code=status.HTTP_204_NO_CONTENT
)
async def delete_task_from_workflow(
    _id: int,
    task_id: int,
    user: User = Depends(current_active_user),
    db: AsyncSession = Depends(get_db),
):
    raise NotImplementedError


@router.patch("{_id}", response_model=WorkflowRead)
async def patch_workflow(
    _id: int,
    patch: WorkflowUpdate,
    user: User = Depends(current_active_user),
    db: AsyncSession = Depends(get_db),
):
    raise NotImplementedError
