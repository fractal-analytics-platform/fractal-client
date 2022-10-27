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


async def test_task_create(db, client, MockCurrentUser):
    """
    GIVEN a task package with a valid __FRACTAL_MANIFEST__
    WHEN it is passed to `POST /task/pypi/`
    THEN
        * the task collection is started in the background
        * at its end the tasks are correctly registered in the database
    """
    raise NotImplementedError
