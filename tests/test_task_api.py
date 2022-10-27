from devtools import debug

from fractal_server.app.models import TaskCreate


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
        command="cmd",
        source="/source",
        input_type="Any",
        output_type="Any",
    )
    async with MockCurrentUser(persist=True):
        res = await client.post("api/v1/task/", json=task.dict())
        debug(res.json())
        data = res.json()
        for key, item in task.dict().items():
            assert data[key] == item
