import pytest
from devtools import debug


async def test_task_collect(clear_db, testserver, register_user, invoke):
    PACKAGE_NAME = "devtools"
    res = await invoke(f"task collect {PACKAGE_NAME}")
    debug(res)
    res.show()


async def test_task_list(clear_db, testserver, register_user, invoke):
    raise NotImplementedError
    # TODO add tasks
    res = await invoke("task list")
    res.show()
    assert res.retcode == 0
    assert len(res.data) == 2


@pytest.mark.xfail
async def test_task_apply(clear_db, testserver, register_user, invoke):
    # TODO:
    # Create project
    # add input resource
    # check that fractal_tasks_core is present (should be optional)
    # try to apply the dummy or another taks
    assert False


async def test_edit_task(
    clear_db, testserver, register_user, invoke, clear_task_cache
):
    raise NotImplementedError
    # TODO create task
    task_name = "FOO"

    res = await invoke(f"task edit {task_name} --name 'new task name'")
    res.show()
    assert res.retcode == 0
    assert res.data["name"] == "new task name"
