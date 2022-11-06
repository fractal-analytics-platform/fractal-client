import asyncio

import pytest
from devtools import debug


async def test_task_collect(
    clear_db, testserver, register_user, invoke, testdata_path
):
    PACKAGE_NAME = testdata_path / "fractal_tasks_dummy-0.1.0-py3-none-any.whl"

    res0 = await invoke(f"task collect {PACKAGE_NAME}")
    debug(res0)
    res0.show()
    venv_path = res0.data["venv_path"]

    res1 = await invoke(f"task check-collection {venv_path}")
    debug(res1)
    res1.show()
    assert res1.data["status"] == "pending"

    await asyncio.sleep(3)

    res2 = await invoke(f"task check-collection {venv_path}")
    debug(res2)
    res2.show()
    assert res2.data["status"] == "OK"


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
