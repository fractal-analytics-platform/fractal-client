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
