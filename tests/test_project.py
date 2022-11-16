from devtools import debug


async def test_project_create(register_user, invoke):
    PROJECT_NAME = "project_name"
    PROJECT_PATH = "project_path"
    res = await invoke(f"project new {PROJECT_NAME} {PROJECT_PATH}")
    debug(res)
    assert res.data["name"] == PROJECT_NAME
    assert res.data["project_dir"] == PROJECT_PATH


async def test_project_create_batch(register_user, invoke):
    res = await invoke("--batch project new project_name project_path")
    debug(res)
    debug(res.data)
    project_id, dataset_id = map(int, res.data.split())
    assert project_id == 1
    assert dataset_id == 1


async def test_project_list(register_user, invoke):
    res = await invoke("project list")
    debug(res)
    debug(vars(res.data))
    assert len(res.data.rows) == 0

    res.show()

    await invoke("--batch project new prj0 prj_path0")
    await invoke("--batch project new prj1 prj_path1")

    res = await invoke("project list")
    debug(res)
    debug(vars(res.data))
    res.show()
    assert len(res.data.rows) == 2


async def test_add_dataset(register_user, invoke):
    DATASET_NAME = "new_ds_name"

    res = await invoke("--batch project new prj0 prj_path0")
    assert res.retcode == 0
    debug(res.data)
    project_id, dataset_id = map(int, res.data.split())

    res = await invoke(f"project add-dataset {project_id} {DATASET_NAME}")
    assert res.retcode == 0
    res.show()
    assert res.data["name"] == DATASET_NAME
