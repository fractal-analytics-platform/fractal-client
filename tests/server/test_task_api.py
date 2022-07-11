from devtools import debug
from sqlmodel import select

from fractal.server.app.models.task import Task
from fractal.server.tasks import collect_tasks


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
        assert data["inserted"] == 1
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
        assert data["updated"] == 1

    res = await db.execute(select(Task))
    n_tasks = len(res.scalars().all())  # FIXME: run query server side!
    assert n_tasks == n_target
