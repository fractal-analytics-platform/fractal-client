from typing import List

from fastapi import APIRouter
from fastapi import BackgroundTasks
from fastapi import Depends
from fastapi import HTTPException
from fastapi import status
from sqlmodel import select

from ...db import AsyncSession
from ...db import get_db
from ...models.project import Dataset
from ...models.project import Project
from ...models.project import ProjectCreate
from ...models.project import ProjectRead
from ...models.project import Resource
from ...models.project import ResourceCreate
from ...models.project import ResourceRead
from ...models.task import Task
from ...runner import submit_workflow
from ...runner import validate_workflow_compatibility
from ...runner import auto_output_dataset
from ...security import current_active_user
from ...security import User
from fractal.common.models import ApplyWorkflow

router = APIRouter()


@router.get("/", response_model=List[ProjectRead])
async def get_list_project(
    user: User = Depends(current_active_user),
    db: AsyncSession = Depends(get_db),
):
    stm = select(Project).where(Project.user_owner_id == user.id)
    res = await db.execute(stm)
    project_list = res.scalars().all()
    return project_list


@router.get("/{project_id}", response_model=ProjectRead)
async def get_project(
    project_id: int,
    user: User = Depends(current_active_user),
    db: AsyncSession = Depends(get_db),
):
    project = await db.get(Project, project_id)
    if not project:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND, detail="Project not found"
        )
    if project.user_owner_id != user.id:
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="Not allowed on project",
        )
    return project


@router.post("/", response_model=ProjectRead, status_code=201)
async def create_project(
    project: ProjectCreate,
    user: User = Depends(current_active_user),
    db: AsyncSession = Depends(get_db),
):
    """
    Create new poject
    """
    db_project = Project.from_orm(project)
    db_project.dataset_list.append(Dataset(name=project.default_dataset_name))

    db_project.user_owner_id = user.id
    db.add(db_project)
    await db.commit()
    await db.refresh(db_project)
    return db_project


@router.post("/{project_id}")
async def add_dataset(
    project_id: int,
):
    """
    Add new dataset to current project
    """
    raise NotImplementedError


@router.post(
    "/{project_id}/{dataset_id}",
    response_model=ResourceRead,
    status_code=status.HTTP_201_CREATED,
)
async def add_resource(
    project_id: int,
    dataset_id: int,
    resource: ResourceCreate,
    user: User = Depends(current_active_user),
    db: AsyncSession = Depends(get_db),
):
    """
    Add resource to an existing dataset
    """
    stm = (
        select(Project, Dataset)
        .join(Dataset)
        .where(Project.user_owner_id == user.id)
    )
    project, dataset = (await db.execute(stm)).one()
    db_resource = Resource(dataset_id=dataset.id, **resource.dict())
    db.add(db_resource)
    await db.commit()
    await db.refresh(db_resource)
    return db_resource


@router.patch("/{project_id}/{dataset_id}")
async def modify_dataset(
    project_id: int,
    dataset_id: int,
):
    raise NotImplementedError


@router.post(
    "/apply/",
    status_code=status.HTTP_202_ACCEPTED,
)
async def apply_workflow(
    apply_workflow: ApplyWorkflow,
    background_tasks: BackgroundTasks,
    user: User = Depends(current_active_user),
    db: AsyncSession = Depends(get_db),
):
    stm = (
        select(Project, Dataset)
        .join(Dataset)
        .where(Project.user_owner_id == user.id)
        .where(Project.id == apply_workflow.project_id)
        .where(Dataset.id == apply_workflow.input_dataset_id)
    )
    project, input_dataset = (await db.execute(stm)).one()

    # TODO check that user is allowed to use this task

    workflow = await db.get(Task, apply_workflow.workflow_id)

    if apply_workflow.output_dataset_id:
        stm = (
            select(Dataset)
            .where(Dataset.project_id == project.id)
            .where(Dataset.id == apply_workflow.output_dataset_id)
        )
        output_dataset = (await db.execute(stm)).one()
    else:
        output_dataset = await auto_output_dataset(
            project=project,
            input_dataset=input_dataset,
            workflow=workflow
        )

    try:
        validate_workflow_compatibility(
            workflow=workflow,
            input_dataset=input_dataset,
            output_dataset=output_dataset
        )
    except TypeError as e:
        raise HTTPException(
            status_code=status.HTTP_422_UNPROCESSABLE_ENTITY, detail=str(e)
        )

    assert input_dataset
    assert output_dataset
    assert workflow

    background_tasks.add_task(
        submit_workflow,
        workflow=workflow,
        input_dataset=input_dataset,
        output_dataset=output_dataset,
    )

    # TODO we should return a job id of some sort
    return dict(status="submitted")
