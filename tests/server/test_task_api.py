import pytest
from devtools import debug
from sqlmodel import select

from fractal.server.app.db import AsyncSession
from fractal.server.app.models.task import Task
from fractal.server.tasks import collect_tasks


@pytest.fixture
async def task_factory():
    """
    Insert task in db
    """

    async def __task_factory(db: AsyncSession, index: int = 0, children=None):
        t = Task(
            name=f"task{index}",
            resource_type="task",
            module=f"task{index}",
            input_type="zarr",
            output_type="zarr",
        )
        if children:
            t.subtask_list.extend(list(children))
        db.add(t)
        await db.commit()
        await db.refresh(t)
        return t

    return __task_factory


async def test_task_relations(db, task_factory):
    t0 = Task(
        name="task0",
        resource_type="task",
        module="task0",
        input_type="zarr",
        output_type="zarr",
    )
    t1 = Task(
        name="task1",
        resource_type="task",
        module="task1",
        input_type="zarr",
        output_type="zarr",
    )
    db.add(t0)
    db.add(t1)
    t0.subtask_list.append(t1)
    await db.commit()
    await db.refresh(t0)
    await db.refresh(t1)
    debug(t0)
    debug(t1)
    assert len(t0.subtask_list) == 1
    assert t0.subtask_list[0] == t1
    assert t0.subtask_list[0].id == t1.id


async def test_collection(db, client, MockCurrentUser):
    """
    GIVEN a running server
    WHEN the `POST task/collect/` endpoint is called
    THEN the table `Task` is updated accordingly
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
        assert data["inserted"] == 2
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
        assert data["updated"] == 2

    res = await db.execute(select(Task))
    n_tasks = len(res.scalars().all())  # FIXME: run query server side!
    assert n_tasks == n_target


async def test_task_get_list(db, client, task_factory, MockCurrentUser):
    t1 = await task_factory(db, index=1)
    await task_factory(db, index=2, children=[t1])

    async with MockCurrentUser(persist=True):
        res = await client.get("api/v1/task/")
        data = res.json()
        assert res.status_code == 200
        debug(data)
        assert len(data) == 2
        assert data[1]["id"] == 2
        assert data[1]["subtask_list"][0]["id"] == 1
