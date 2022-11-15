from devtools import debug  # noqa
from sqlmodel import select

from fractal_server.app.models import TaskRead
from fractal_server.app.models import Workflow
from fractal_server.app.models import WorkflowRead
from fractal_server.app.models import WorkflowTask


async def test_post_workflow(db, client, MockCurrentUser, project_factory):
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


async def test_delete_workflow(
    db, client, MockCurrentUser, project_factory, task_factory
):
    """
    GIVEN a Workflow with two Tasks
    WHEN the endpoint that DELETE a Workflow is called
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


async def test_get_workflow(
    db, client, MockCurrentUser, project_factory, task_factory
):
    """
    GIVEN a Workflow in the db
    WHEN the endpoint to GET a Workflow by its id is called
    THEN the Workflow is returned
    """
    async with MockCurrentUser(persist=True) as user:
        project = await project_factory(user)
        p_id = project.id
        workflow = {
            "name": "My Workflow",
            "project_id": p_id,
            "task_list": [],
        }
        res = await client.post(
            "api/v1/workflow/",
            json=workflow,
        )
        wf_id = res.json()["id"]
        res = await client.get(f"/api/v1/workflow/{wf_id}")

        assert res.status_code == 200
        workflow["id"] = wf_id
        assert res.json() == workflow


async def test_post_newtask(
    db, client, MockCurrentUser, project_factory, task_factory
):
    """
    GIVEN a Workflow with a list of WorkflowTasks
    WHEN the endpoint to POST a new WorkflowTask inside
        the Workflow.task_list is called
    THEN the new WorkflowTask is inserted in Workflow.task_list
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
        assert res.status_code == 201
        wf_id = res.json()["id"]

        workflow = await db.get(Workflow, wf_id)
        t0 = await task_factory()
        t1 = await task_factory()
        await workflow.insert_task(t0.id, db=db)
        await workflow.insert_task(t1.id, db=db)
        await db.refresh(workflow)

        assert len(workflow.task_list) == 2
        assert workflow.task_list[0].task == t0
        assert workflow.task_list[1].task == t1

        await db.refresh(workflow)

        t2 = await task_factory()
        last_task = {"task_id": t2.id, "args": {"a": 0, "b": 1}}

        res = await client.post(
            f"api/v1/workflow/{wf_id}/add-task/",
            json=last_task,
        )
        assert res.status_code == 201

        t0b = await task_factory()
        second_task = {
            "task_id": t0b.id,
            "order": 1,
        }
        res = await client.post(
            f"api/v1/workflow/{wf_id}/add-task/",
            json=second_task,
        )
        assert res.status_code == 201

        workflow = WorkflowRead(**res.json())

        assert len(workflow.task_list) == 4
        assert workflow.task_list[0].task == TaskRead(**t0.dict())
        assert workflow.task_list[1].task == TaskRead(**t0b.dict())
        assert workflow.task_list[2].task == TaskRead(**t1.dict())
        assert workflow.task_list[3].task == TaskRead(**t2.dict())
        assert workflow.task_list[3].args == last_task["args"]


async def test_delete_task(
    db, client, MockCurrentUser, project_factory, task_factory
):
    """
    GIVEN a Workflow with a list of WorkflowTasks
    WHEN the endpoint to DELETE a WorkflowTask in the
        Workflow.task_list is called
    THEN the selected WorkflowTask is properly removed
        from Workflow.task_list
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
        assert res.status_code == 201
        wf_id = res.json()["id"]

        workflow = await db.get(Workflow, wf_id)
        t0 = await task_factory()
        t1 = await task_factory()
        t2 = await task_factory()

        await workflow.insert_task(t0.id, db=db)
        await workflow.insert_task(t1.id, db=db)
        await workflow.insert_task(t2.id, db=db)
        await db.refresh(workflow)

        assert (
            len((await db.execute(select(WorkflowTask))).scalars().all()) == 3
        )
        assert len(workflow.task_list) == 3
        for i, task in enumerate(workflow.task_list):
            assert task.order == i

        res = await client.delete(f"api/v1/workflow/{wf_id}/rm-task/{t1.id}")
        assert res.status_code == 204

        await db.refresh(workflow)
        assert (
            len((await db.execute(select(WorkflowTask))).scalars().all()) == 2
        )
        assert len(workflow.task_list) == 2
        for i, task in enumerate(workflow.task_list):
            assert task.order == i


async def test_get_project_workflows(
    db, client, MockCurrentUser, project_factory, task_factory
):
    """
    GIVEN a Project containing three Workflows
    WHEN the endpoint to GET all the Workflows associated
        to that Project is called
    THEN the list of all its Workflows is returned
    """
    async with MockCurrentUser(persist=True) as user:
        project = await project_factory(user)
        other_project = await project_factory(user)
        workflow1 = {"name": "WF1", "project_id": project.id}
        workflow2 = {"name": "WF2", "project_id": project.id}
        workflow3 = {"name": "WF3", "project_id": other_project.id}
        workflow4 = {"name": "WF4", "project_id": project.id}
        res = await client.post("api/v1/workflow/", json=workflow1)
        res = await client.post("api/v1/workflow/", json=workflow2)
        res = await client.post("api/v1/workflow/", json=workflow3)
        res = await client.post("api/v1/workflow/", json=workflow4)

        res = await client.get(f"api/v1/project/{project.id}/workflows/")

        workflow_list = res.json()
        assert len(workflow_list) == 3
        assert len((await db.execute(select(Workflow))).scalars().all()) == 4


async def test_patch_workflow(
    db, client, MockCurrentUser, project_factory, task_factory
):
    """
    GIVEN a Workflow
    WHEN the endpoint to PATCH a Workflow is called
    THEN the Workflow is updated
    """
    async with MockCurrentUser(persist=True) as user:
        project = await project_factory(user)
        other_project = await project_factory(user)
        workflow = {"name": "WF", "project_id": project.id}
        res = await client.post("api/v1/workflow/", json=workflow)
        wf_id = res.json()["id"]
        assert res.json()["name"] == "WF"

        workflow = await db.get(Workflow, wf_id)
        res = await client.get(f"api/v1/project/{project.id}/workflows/")
        assert len(res.json()) == 1
        res = await client.get(f"api/v1/project/{other_project.id}/workflows/")
        assert len(res.json()) == 0

        patch = {"name": "FW", "project_id": other_project.id}
        res = await client.patch(f"api/v1/workflow/{wf_id}", json=patch)

        await db.refresh(workflow)
        assert workflow.name == "FW"
        res = await client.get(f"api/v1/project/{project.id}/workflows/")
        assert len(res.json()) == 0
        res = await client.get(f"api/v1/project/{other_project.id}/workflows/")
        assert len(res.json()) == 1


async def test_patch_workflow_task(
    db, client, MockCurrentUser, project_factory, task_factory
):
    """
    GIVEN a WorkflowTask
    WHEN the endpoint to PATCH a WorkflowTask is called
    THEN the WorkflowTask is updated
    """
    async with MockCurrentUser(persist=True) as user:
        project = await project_factory(user)
        workflow = {"name": "WF", "project_id": project.id}
        res = await client.post("api/v1/workflow/", json=workflow)
        wf_id = res.json()["id"]

        workflow = await db.get(Workflow, wf_id)
        t0 = await task_factory()
        await workflow.insert_task(t0.id, db=db)
        await db.refresh(workflow)
        assert workflow.task_list[0].args is None

        payload = dict(args={"a": 123, "d": 321}, executor="cpu-low")
        res = await client.patch(
            f"api/v1/workflow/{workflow.id}/"
            f"edit-task/{workflow.task_list[0].id}",
            json=payload,
        )

        patched_workflow_task = res.json()
        assert patched_workflow_task["args"] == payload["args"]
        assert patched_workflow_task["executor"] == payload["executor"]
        assert res.status_code == 200

        payload_up = dict(args={"a": {"c": 43}, "b": 123})
        res = await client.patch(
            f"api/v1/workflow/{workflow.id}/"
            f"edit-task/{workflow.task_list[0].id}",
            json=payload_up,
        )
        patched_workflow_task_up = res.json()
        assert patched_workflow_task_up["args"] == dict(
            a=dict(c=43), b=123, d=321
        )
