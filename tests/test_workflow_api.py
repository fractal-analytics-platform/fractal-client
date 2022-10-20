from sqlmodel import select

from fractal_server.app.models import Workflow
from fractal_server.app.models import WorkflowTask


async def test_workflow_post(db, client, MockCurrentUser, project_factory):
    async with MockCurrentUser(persist=True) as user:
        res = await client.post(
            "api/v1/workflow/",
            json={"name": "My Workflow"},
        )
        assert res.status_code == 422  # no project_id

        project1 = await project_factory(user)
        p1_id = project1.id
        workflow1 = {
            "name": "My Workflow",
            "project_id": p1_id,
        }
        project2 = await project_factory(user)
        p2_id = project2.id
        workflow2 = {
            "name": "My Workflow",
            "project_id": p2_id,
        }

        res = await client.post(
            "api/v1/workflow/",
            json=workflow1,
        )
        assert res.status_code == 201

        res = await client.post(
            "api/v1/workflow/",
            json=workflow1,
        )
        assert res.status_code == 422  # already in use

        res = await client.post(
            "api/v1/workflow/",
            json={
                "name": "My Workflow",
                "project_id": p1_id + p2_id,
            },
        )
        assert res.status_code == 404  # project does not exist

        res = await client.post(
            "api/v1/workflow/",
            json=workflow2,
        )
        assert res.status_code == 201  # same name, different projects

        for id in [p1_id, p2_id]:
            stm = select(Workflow).where(Workflow.project_id == id)
            _workflow = await db.execute(stm)
            db_workflow = _workflow.scalars().one()

            assert db_workflow.name == "My Workflow"
            assert db_workflow.project_id == id


async def test_workflow_delete(
    db, client, MockCurrentUser, project_factory, task_factory
):
    """
    GIVEN a Workflow with two Tasks
    WHEN the delete endpoint is called
    THEN the Workflow and its associated WorkflowTasks
        are removed from the db
    """
    async with MockCurrentUser(persist=True) as user:
        project = await project_factory(user)
        p_id = project.id
        workflow = {
            "name": "My Workflow",
            "project_id": p_id,
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
        assert len((await db.execute(wf_tasks)).scalars().all()) == 2

        res = await client.delete(f"api/v1/workflow/{wf_id}")
        assert res.status_code == 204

        assert (await db.get(Workflow, wf_id)) is None
        assert len((await db.execute(wf_tasks)).scalars().all()) == 0
