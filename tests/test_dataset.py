async def test_add_resource(clear_db, testserver, register_user, invoke):
    res = await invoke("project new prj0 prj_path0")
    project_id = res.data["id"]
    dataset_id = res.data["dataset_list"][0]["id"]
    assert res.data["dataset_list"][0]["resource_list"] == []

    PATH = "/new/ds/path"
    res = await invoke(
        f"dataset add-resource {project_id} {dataset_id} {PATH}"
    )
    res.show()
    assert res.retcode == 0
    assert res.data["path"] == PATH
    assert res.data["dataset_id"] == dataset_id


async def test_edit_dataset(clear_db, testserver, register_user, invoke):
    res = await invoke("project new prj0 prj_path0")
    project_id = res.data["id"]
    dataset_id = res.data["dataset_list"][0]["id"]

    res = await invoke(f"dataset edit {project_id} {dataset_id}")
