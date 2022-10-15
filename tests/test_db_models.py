import pytest
from devtools import debug
from sqlalchemy.exc import IntegrityError
from sqlmodel import select

from fractal_server.app.models import Project
from fractal_server.app.models.workflow import Workflow
from fractal_server.app.models.workflow import WorkflowTask


async def test_project_name_not_unique(MockCurrentUser, db, project_factory):
    """
    GIVEN the fractal_server database
    WHEN I create two projects with the same name and same user
    THEN no exception is raised
    """
    PROJ_NAME = "project name"
    async with MockCurrentUser(persist=True) as user:
        p0 = await project_factory(user, name=PROJ_NAME)
        p1 = await project_factory(user, name=PROJ_NAME)

    stm = select(Project).where(Project.name == PROJ_NAME)
    res = await db.execute(stm)
    project_list = res.scalars().all()
    assert len(project_list) == 2
    assert p0 in project_list
    assert p1 in project_list


async def test_project_workflow_relation(db, project_factory, MockCurrentUser):
    """
    GIVEN no preconditions
    WHEN I create a workflow
    THEN I must associate the workflow to a project or error
    """

    # Fail if no project_id
    with pytest.raises(IntegrityError):
        wf = Workflow(name="my wfl")
        db.add(wf)
        await db.commit()
    await db.rollback()

    async with MockCurrentUser(persist=True) as user:
        project = await project_factory(user)
        wf = Workflow(name="my wfl", project_id=project.id)
        db.add(wf)
        await db.commit()


async def test_task_workflow_association(
    db, project_factory, MockCurrentUser, task_factory
):
    async with MockCurrentUser(persist=True) as user:
        project = await project_factory(user)
        t0 = await task_factory()
        t1 = await task_factory()

        wf = Workflow(name="my wfl", project_id=project.id)
        await wf.insert_task(t0, db=db)

        db.add(wf)
        await db.commit()
        await db.refresh(wf)

        debug(wf)
        # check workflow
        assert len(wf.task_list) == 1
        assert wf.task_list[0].id == t0.id
        # check association table
        stm = (
            select(WorkflowTask)
            .where(WorkflowTask.workflow_id == wf.id)
            .where(WorkflowTask.task_id == t0.id)
        )
        res = await db.execute(stm)
        link = res.scalars().one()
        debug(link)
        assert link.task_id == t0.id

        # Insert at position 0
        await wf.insert_task(t1, order=0, db=db)
        db.add(wf)
        await db.commit()
        await db.refresh(wf)

        stm = (
            select(WorkflowTask)
            .where(WorkflowTask.workflow_id == wf.id)
            .where(WorkflowTask.task_id == t1.id)
        )
        res = await db.execute(stm)
        link = res.scalars().one()
        debug(link)
        assert link.order == 0
        assert link.task_id == t1.id
