import asyncio
import logging
from pathlib import Path
from typing import List

from fastapi import APIRouter
from fastapi import BackgroundTasks
from fastapi import Depends
from fastapi import HTTPException
from fastapi import Response
from fastapi import status
from pydantic import UUID4
from sqlalchemy.exc import IntegrityError
from sqlmodel import select

from ...db import AsyncSession
from ...db import DBSyncSession
from ...db import get_db
from ...db import get_sync_db
from ...models import ApplyWorkflow
from ...models import ApplyWorkflowCreate
from ...models import ApplyWorkflowRead
from ...models import Dataset
from ...models import DatasetCreate
from ...models import DatasetRead
from ...models import DatasetUpdate
from ...models import LinkUserProject
from ...models import Project
from ...models import ProjectCreate
from ...models import ProjectRead
from ...models import Resource
from ...models import ResourceCreate
from ...models import ResourceRead
from ...models import ResourceUpdate
from ...models.task import Task
from ...runner import auto_output_dataset
from ...runner import submit_workflow
from ...runner import validate_workflow_compatibility
from ...security import current_active_user
from ...security import User


router = APIRouter()


async def get_project_check_owner(
    *,
    project_id: int,
    user_id: UUID4,
    db: AsyncSession = Depends(get_db),
) -> Project:
    """
    Check that user is a member of project and return

    Raise 403_FORBIDDEN if the user is not a member
    Raise 404_NOT_FOUND if the project does not exist
    """
    project, link_user_project = await asyncio.gather(
        db.get(Project, project_id),
        db.get(LinkUserProject, (project_id, user_id)),
    )
    if not project:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND, detail="Project not found"
        )
    if not link_user_project:
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail=f"Not allowed on project {project_id}",
        )
    return project


@router.get("/", response_model=List[ProjectRead])
async def get_list_project(
    user: User = Depends(current_active_user),
    db: AsyncSession = Depends(get_db),
) -> List[Project]:
    """
    Return list of projects user is member of
    """
    stm = (
        select(Project)
        .join(LinkUserProject)
        .where(LinkUserProject.user_id == user.id)
    )
    res = await db.execute(stm)
    project_list = res.scalars().all()
    return project_list


@router.get("/{project_id}", response_model=ProjectRead)
async def get_project(
    project_id: int,
    user: User = Depends(current_active_user),
    db: AsyncSession = Depends(get_db),
):
    project = await get_project_check_owner(
        project_id=project_id, user_id=user.id, db=db
    )
    return project


@router.delete("/{project_id}", status_code=204)
async def delete_project(
    project_id: int,
    user: User = Depends(current_active_user),
    db: AsyncSession = Depends(get_db),
):
    project = await get_project_check_owner(
        project_id=project_id, user_id=user.id, db=db
    )
    await db.delete(project)
    await db.commit()
    return Response(status_code=status.HTTP_204_NO_CONTENT)


@router.post("/", response_model=ProjectRead, status_code=201)
async def create_project(
    project: ProjectCreate,
    user: User = Depends(current_active_user),
    db: AsyncSession = Depends(get_db),
):
    """
    Create new poject
    """

    # Check that there is no project with the same user and name
    stm = (
        select(Project)
        .join(LinkUserProject)
        .where(Project.name == project.name)
        .where(LinkUserProject.user_id == user.id)
    )
    res = await db.execute(stm)
    if res.scalars().all():
        raise HTTPException(
            status_code=status.HTTP_422_UNPROCESSABLE_ENTITY,
            detail=f"Project name ({project.name}) already in use",
        )

    db_project = Project.from_orm(project)
    db_project.dataset_list.append(Dataset(name=project.default_dataset_name))
    db_project.user_member_list.append(user)
    try:
        db.add(db_project)
        await db.commit()
        await db.refresh(db_project)
    except IntegrityError as e:
        await db.rollback()
        logging.error(str(e))
        raise HTTPException(
            status_code=status.HTTP_422_UNPROCESSABLE_ENTITY,
            detail=str(e),
        )

    return db_project


@router.post(
    "/apply/",
    status_code=status.HTTP_202_ACCEPTED,
    response_model=ApplyWorkflowRead,
)
async def apply_workflow(
    apply_workflow: ApplyWorkflowCreate,
    background_tasks: BackgroundTasks,
    user: User = Depends(current_active_user),
    db: AsyncSession = Depends(get_db),
    db_sync: DBSyncSession = Depends(get_sync_db),
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

    workflow = db_sync.get(Task, apply_workflow.workflow_id)
    if not workflow:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Workflow {apply_workflow.workflow_id} not found",
        )

    if apply_workflow.output_dataset_id:
        stm = (
            select(Dataset)
            .where(Dataset.project_id == project.id)
            .where(Dataset.id == apply_workflow.output_dataset_id)
        )
        output_dataset = (await db.execute(stm)).scalars().one()
    else:
        output_dataset = await auto_output_dataset(
            project=project, input_dataset=input_dataset, workflow=workflow
        )

    try:
        validate_workflow_compatibility(
            workflow=workflow,
            input_dataset=input_dataset,
            output_dataset=output_dataset,
        )
    except TypeError as e:
        raise HTTPException(
            status_code=status.HTTP_422_UNPROCESSABLE_ENTITY, detail=str(e)
        )

    if not input_dataset:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Dataset {apply_workflow.dataset_id} not found",
        )
    if not output_dataset:
        raise HTTPException(
            status_code=status.HTTP_422_UNPROCESSABLE_ENTITY,
            detail="Could not determine output dataset.",
        )

    job = ApplyWorkflow.from_orm(apply_workflow)
    db.add(job)
    await db.commit()
    await db.refresh(job)

    background_tasks.add_task(
        submit_workflow,
        job=job,
        user=user,
        db=db,
    )

    return job


@router.post(
    "/{project_id}/",
    response_model=DatasetRead,
    status_code=status.HTTP_201_CREATED,
)
async def add_dataset(
    project_id: int,
    dataset: DatasetCreate,
    user: User = Depends(current_active_user),
    db: AsyncSession = Depends(get_db),
):
    """
    Add new dataset to current project
    """
    project = await get_project_check_owner(
        project_id=project_id, user_id=user.id, db=db
    )
    dataset.project_id = project.id
    db_dataset = Dataset.from_orm(dataset)
    db.add(db_dataset)
    await db.commit()
    await db.refresh(db_dataset)
    return db_dataset


@router.delete("/{project_id}/{dataset_id}", status_code=204)
async def delete_dataset(
    project_id: int,
    dataset_id: int,
    user: User = Depends(current_active_user),
    db: AsyncSession = Depends(get_db),
):
    await get_project_check_owner(
        project_id=project_id, user_id=user.id, db=db
    )
    stm = (
        select(Dataset)
        .join(Project)
        .where(Project.id == project_id)
        .where(Dataset.id == dataset_id)
    )
    res = await db.execute(stm)
    dataset = res.scalar()
    await db.delete(dataset)
    await db.commit()
    return Response(status_code=status.HTTP_204_NO_CONTENT)


@router.get(
    "/{project_id}/{dataset_id}",
    response_model=List[ResourceRead],
)
async def get_resource(
    project_id: int,
    dataset_id: int,
    user: User = Depends(current_active_user),
    db: AsyncSession = Depends(get_db),
):
    """
    Get resources from a dataset
    """
    await get_project_check_owner(
        project_id=project_id, user_id=user.id, db=db
    )
    stm = select(Resource).where(Resource.dataset_id == dataset_id)
    res = await db.execute(stm)
    resource_list = res.scalars().all()
    return resource_list


@router.delete("/{project_id}/{dataset_id}/{resource_id}", status_code=204)
async def delete_resource(
    project_id: int,
    dataset_id: int,
    resource_id: int,
    user: User = Depends(current_active_user),
    db: AsyncSession = Depends(get_db),
):
    project = await get_project_check_owner(
        project_id=project_id, user_id=user.id, db=db
    )
    resource = await db.get(Resource, resource_id)
    if not resource or resource.dataset_id not in (
        ds.id for ds in project.dataset_list
    ):
        raise HTTPException(
            status_code=status.HTTP_422_UNPROCESSABLE_ENTITY,
            detail="Resource does not exist or does not belong to project",
        )
    await db.delete(resource)
    await db.commit()
    return Response(status_code=status.HTTP_204_NO_CONTENT)


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

    # Check that path is absolute, which is needed for when the server submits
    # tasks as a different user

    if not Path(resource.path).is_absolute():
        raise HTTPException(
            status_code=status.HTTP_422_UNPROCESSABLE_ENTITY,
            detail=f"Path `{resource.path}` is not absolute.",
        )

    project = await get_project_check_owner(
        project_id=project_id, user_id=user.id, db=db
    )
    dataset = await db.get(Dataset, dataset_id)
    if dataset not in project.dataset_list:
        raise HTTPException(
            status_code=status.HTTP_422_UNPROCESSABLE_ENTITY,
            detail=f"Dataset {dataset_id} is not part of project {project_id}",
        )

    db_resource = Resource(dataset_id=dataset.id, **resource.dict())
    db.add(db_resource)
    await db.commit()
    await db.refresh(db_resource)
    return db_resource


@router.patch(
    "/{project_id}/{dataset_id}/{resource_id}", response_model=ResourceRead
)
async def edit_resource(
    project_id: int,
    dataset_id: int,
    resource_id: int,
    resource_update: ResourceUpdate,
    user: User = Depends(current_active_user),
    db: AsyncSession = Depends(get_db),
):
    project = await get_project_check_owner(
        project_id=project_id, user_id=user.id, db=db
    )
    dataset = await db.get(Dataset, dataset_id)
    orig_resource = await db.get(Resource, resource_id)

    if dataset not in project.dataset_list:
        raise HTTPException(
            status_code=status.HTTP_422_UNPROCESSABLE_ENTITY,
            detail=f"Dataset {dataset_id} is not part of project {project_id}",
        )
    if orig_resource not in dataset.resource_list:
        raise HTTPException(
            status_code=status.HTTP_422_UNPROCESSABLE_ENTITY,
            detail=(
                f"Resource {resource_id} is not part of "
                f"dataset {dataset_id}"
            ),
        )

    for key, value in resource_update.dict(exclude_unset=True).items():
        setattr(orig_resource, key, value)
    await db.commit()
    await db.refresh(orig_resource)
    return orig_resource


@router.patch("/{project_id}/{dataset_id}", response_model=DatasetRead)
async def patch_dataset(
    project_id: int,
    dataset_id: int,
    dataset_update: DatasetUpdate,
    user: User = Depends(current_active_user),
    db: AsyncSession = Depends(get_db),
):
    project = await get_project_check_owner(
        project_id=project_id, user_id=user.id, db=db
    )
    db_dataset = await db.get(Dataset, dataset_id)
    if db_dataset not in project.dataset_list:
        raise HTTPException(
            status_code=status.HTTP_422_UNPROCESSABLE_ENTITY,
            detail=f"Dataset {dataset_id} is not part of project {project_id}",
        )

    for key, value in dataset_update.dict(exclude_unset=True).items():
        setattr(db_dataset, key, value)

    await db.commit()
    await db.refresh(db_dataset)
    return db_dataset
