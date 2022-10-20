from datetime import timezone

from devtools import debug
from sqlmodel import select

from fractal_server.app.models import Workflow
from fractal_server.app.models import WorkflowTask
from fractal_server.app.models.models_utils import get_timestamp


def test_timestamp():
    """
    GIVEN a function that provides a timestamp
    WHEN called
    THEN the timestamp is timezone aware and the timezone is set to UTC
    """
    ts = get_timestamp()
    debug(ts)
    assert ts
    assert ts.tzinfo is timezone.utc


async def test_cascade_delete_workflow(
    db, client, MockCurrentUser, project_factory, task_factory
):
    """
    GIVEN a Workflow
    WHEN the Workflow is deleted via API
    THEN all the related WorkflowTask are deleted
    """
    async with MockCurrentUser(persist=True) as user:
        project = await project_factory(user)
        workflow = {
            "name": "My Workflow",
            "project_id": project.id,
        }
        res = await client.post(
            "api/v1/workflow/",
            json=workflow,
        )
        wf_id = res.json()["id"]

        wf = await db.get(Workflow, wf_id)
        t0 = await task_factory()
        t1 = await task_factory()
        await wf.insert_task(t0, db=db)
        await wf.insert_task(t1, db=db)

        wf_tasks = select(WorkflowTask).where(
            WorkflowTask.workflow_id == wf_id
        )
        assert len(db.execute(wf_tasks).scalars().all()) == 2

        await client.delete(f"/api/v1/workflow/{wf_id}")

        assert len(db.execute(wf_tasks).scalars().all()) == 0


async def test_cascade_delete_project(
    db, client, MockCurrentUser, project_factory, task_factory
):
    """
    GIVEN a Project
    WHEN the Project is deleted via API
    THEN all the related Workflows are deleted
    """
    async with MockCurrentUser(persist=True) as user:
        project = await project_factory(user)
        workflow1 = {
            "name": "My First Workflow",
            "project_id": project.id,
        }
        await client.post(
            "api/v1/workflow/",
            json=workflow1,
        )
        workflow2 = {
            "name": "My Second Workflow",
            "project_id": project.id,
        }
        await client.post(
            "api/v1/workflow/",
            json=workflow2,
        )

        workflows = select(Workflow).where(Workflow.project_id == project.id)
        assert len(db.execute(workflows).scalars().all()) == 2

        await client.delete(f"/api/v1/project/{project.id}")

        assert len(db.execute(workflows).scalars().all()) == 0
