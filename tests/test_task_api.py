from devtools import debug


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


async def test_collection_api(client, dummy_task_package, MockCurrentUser):
    """
    GIVEN a package in a format that `pip` understands
    WHEN the api to collect tasks from that package is called
    THEN an environment is created, the package is installed and the task
         collected
    """
    PREFIX = "/api/v1/task"

    task_collection = dict(package=dummy_task_package.as_posix())

    async with MockCurrentUser(persist=True):
        res = await client.post(f"{PREFIX}/pip/", json=task_collection)
        debug(res.json())
        assert res.status_code == 201

    task_list = res.json()
    task_names = (t["name"] for t in task_list)
    assert len(task_list) == 2
    assert "dummy" in task_names
    assert "dummy parallel" in task_names
