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


async def test_task_apply(clear_db, testserver, register_user, invoke):
    from devtools import debug

    PROJECT_NAME = "project_name"
    PROJECT_PATH = "project_path"
    res = await invoke(f"project new {PROJECT_NAME} {PROJECT_PATH}")
    project_id = res.data["id"]
    dataset_id = res.data["dataset_list"][0]["id"]
    res = await invoke(f"dataset edit {project_id} {dataset_id} --type image")
    res = await invoke("task new mytask task image zarr mypackage.subpkg:foo")
    res = await invoke("task list")
    res.show()
    res = await invoke("task apply project_name default default mytask ")
    debug(res)
