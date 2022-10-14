import pytest
from devtools import debug
from sqlmodel import select

from fractal_server.app.models import Task
from fractal_server.app.models import TaskCreate
from fractal_server.tasks import collect_tasks

try:
    from fractal_tasks_core import __FRACTAL_MANIFEST__

    HAS_TASK_CORE = True
    N_CORE_TASKS = len(__FRACTAL_MANIFEST__)
except ImportError:
    HAS_TASK_CORE = False


@pytest.mark.skipif(
    not HAS_TASK_CORE, reason="fractal-tasks-core not installed"
)
async def test_collection(db, client, MockCurrentUser):
    """
    GIVEN a running server
    WHEN the `POST task/collect/` endpoint is called
    THEN the table `Task` is updated accordingly, collecting the available
         tasks
    """
    res = await db.execute(select(Task))
    n_tasks = len(res.scalars().all())  # FIXME: run query server side!
    assert n_tasks == 0

    n_target = len(list(collect_tasks()))

    async with MockCurrentUser(persist=True):
        res = await client.post("api/v1/task/collect/")
        assert res.status_code == 201
        data = res.json()
        debug(data)
        assert data["inserted"] == N_CORE_TASKS
        assert data["updated"] == 0

    res = await db.execute(select(Task))
    n_tasks = len(res.scalars().all())  # FIXME: run query server side!
    assert n_tasks == n_target

    # Check for idempotency
    async with MockCurrentUser(persist=True):
        res = await client.post("api/v1/task/collect/")
        data = res.json()
        assert res.status_code == 201
        assert data["inserted"] == 0
        assert data["updated"] == N_CORE_TASKS

    res = await db.execute(select(Task))
    n_tasks = len(res.scalars().all())  # FIXME: run query server side!
    assert n_tasks == n_target


async def test_task_get_list(db, client, task_factory, MockCurrentUser):
    t0 = await task_factory(name="task0")
    t1 = await task_factory(name="task1")
    t2 = await task_factory(index=2, subtask_list=[t0, t1])

    async with MockCurrentUser(persist=True):
        res = await client.get("api/v1/task/")
        data = res.json()
        assert res.status_code == 200
        debug(data)
        assert len(data) == 3
        assert data[2]["id"] == t2.id


async def test_task_create(db, client, MockCurrentUser):
    """
    GIVEN a CreateTask object
    WHEN it is fed to the `POST task` endpoint
    THEN a new task is correctly created
    """
    task = TaskCreate(
        name="mytask",
        resource_type="workflow",
        input_type="Any",
        output_type="Any",
    )
    async with MockCurrentUser(persist=True):
        res = await client.post("api/v1/task/", json=task.dict())
        debug(res.json())
        data = res.json()
        for key, item in task.dict().items():
            assert data[key] == item
