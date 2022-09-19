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
