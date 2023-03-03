import pytest
from devtools import debug


async def test_project_create(register_user, invoke, tmp_path):
    PROJECT_NAME = "project_name"
    PROJECT_PATH = str(tmp_path)
    res = await invoke(f"project new {PROJECT_NAME} {PROJECT_PATH}")
    debug(res)
    assert res.data["name"] == PROJECT_NAME
    assert res.data["project_dir"] == PROJECT_PATH


async def test_project_delete(register_user, invoke, tmp_path):

    # Create project
    project_dir = str(tmp_path)
    res = await invoke(f"project new MyProj1 {project_dir}")
    res.show()
    project_id_1 = res.data["id"]

    # Show project
    res = await invoke(f"project show {project_id_1}")
    res.show()

    # Delete project
    res = await invoke(f"project delete {project_id_1}")
    assert res.retcode == 0

    # Try to show deleted project, and fail
    with pytest.raises(SystemExit):
        res = await invoke(f"project show {project_id_1}")


async def test_project_create_batch(register_user, invoke, tmp_path):
    project_dir = str(tmp_path)
    res = await invoke(f"--batch project new MyProj1 {project_dir}")
    debug(res)
    debug(res.data)
    project_id, dataset_id = map(int, res.data.split())
    assert project_id == 1
    assert dataset_id == 1


async def test_project_list(register_user, invoke, tmp_path):
    res = await invoke("project list")
    debug(res)
    debug(vars(res.data))
    assert len(res.data.rows) == 0

    res.show()

    project_dir = str(tmp_path)
    res = await invoke(f"--batch project new proj0 {project_dir}")
    res = await invoke(f"--batch project new proj1 {project_dir}")

    res = await invoke("project list")
    debug(res)
    debug(vars(res.data))
    res.show()
    assert len(res.data.rows) == 2


async def test_add_dataset(register_user, invoke, tmp_path):
    DATASET_NAME = "new_ds_name"

    project_dir = str(tmp_path)
    res = await invoke(f"--batch project new proj0 {project_dir}")
    assert res.retcode == 0
    debug(res.data)
    project_id, dataset_id = map(int, res.data.split())

    res = await invoke(f"project add-dataset {project_id} {DATASET_NAME}")
    assert res.retcode == 0
    res.show()
    assert res.data["name"] == DATASET_NAME
