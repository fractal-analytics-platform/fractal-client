from devtools import debug  # noqa
from sqlmodel import select

from fractal_server.app.models import Workflow


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

        res = await client.delete(f"api/v1/workflow/{wf_id}")
        assert res.status_code == 204

        # TODO add tasks with API and test cascade delete

        assert not await db.get(Workflow, wf_id)


async def test_workflow_get(
    db, client, MockCurrentUser, project_factory, task_factory
):
    """
    GIVEN a Workflow
    WHEN the get endpoint is called
    THEN the Workflow is returned
    """
    async with MockCurrentUser(persist=True) as user:
        project = await project_factory(user)
        p_id = project.id
        WF_ID = 1
        workflow = {
            "id": WF_ID,
            "name": "My Workflow",
            "project_id": p_id,
            "task_list": [],
        }
        res = await client.post(
            "api/v1/workflow/",
            json=workflow,
        )

        res = await client.get(f"api/v1/workflow/{WF_ID}")
        assert res.status_code == 200
        assert res.json() == workflow
