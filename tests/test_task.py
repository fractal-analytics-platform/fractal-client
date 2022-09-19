import pytest


async def test_task_new(clear_db, testserver, register_user, invoke):
    res = await invoke("task new mytask task image zarr mypackage.subpkg:foo")
    res.show()
    assert res.retcode == 0
    assert res.data["name"] == "mytask"
    assert res.data["module"] == "mypackage.subpkg:foo"


async def test_task_list(clear_db, testserver, register_user, invoke):
    res = await invoke("task new mytask0 task image zarr mypackage.subpkg:foo")
    res = await invoke("task new mytask1 task image zarr mypackage.subpkg:foo")
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


async def test_edit_task(clear_db, testserver, register_user, invoke):
    res = await invoke("task new mytask0 task image zarr mypackage.subpkg:foo")
    task_id = res.data["id"]

    res = await invoke(f"task edit {task_id} --name 'new task name'")
    res.show()
    assert res.retcode == 0
    assert res.data["name"] == "new task name"


async def test_add_subtask(clear_db, testserver, register_user, invoke):
    res = await invoke("task new parent task image zarr mypackage.subpkg:foo")
    parent_task_id = res.data["id"]
    res = await invoke("task new child task image zarr mypackage.subpkg:foo")
    child_task_id = res.data["id"]

    res = await invoke(f"task add-subtask {parent_task_id} {child_task_id}")
    res.show()

    assert res.retcode == 0
    assert res.data["subtask_list"][0]["subtask_id"] == child_task_id