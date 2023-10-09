import pytest
from devtools import debug


async def test_project_create(register_user, invoke):
    PROJECT_NAME = "project_name"
    res = await invoke(f"project new {PROJECT_NAME}")
    debug(res)
    assert res.data["name"] == PROJECT_NAME


async def test_project_delete(register_user, invoke):

    # Create project
    res = await invoke("project new MyProj1")
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


async def test_project_create_batch(register_user, invoke):
    res = await invoke("--batch project new MyProj1")
    debug(res)
    debug(res.data)
    project_id = int(res.data)
    assert project_id == 1


async def test_project_list(register_user, invoke):
    res = await invoke("project list")
    debug(res)
    debug(vars(res.data))
    assert len(res.data.rows) == 0

    res.show()

    res = await invoke("--batch project new proj0")
    res = await invoke("--batch project new proj1")

    res = await invoke("project list")
    debug(res)
    debug(vars(res.data))
    res.show()
    assert len(res.data.rows) == 2


async def test_add_dataset(register_user, invoke):
    DATASET_NAME = "new_ds_name"

    res = await invoke("--batch project new proj0")
    assert res.retcode == 0
    debug(res.data)
    project_id = int(res.data)

    res = await invoke(f"project add-dataset {project_id} {DATASET_NAME}")
    assert res.retcode == 0
    res.show()
    assert res.data["name"] == DATASET_NAME


@pytest.mark.parametrize("new_name", ["new_name", None])
@pytest.mark.parametrize("read_only", [True, False, None])
async def test_edit_project(
    register_user,
    invoke,
    new_name,
    read_only,
    tmp_path,
):
    name = "name"
    res = await invoke(f"project new {name}")
    project = res.data
    project_id = project["id"]

    cmd = f"project edit {project_id}"
    if new_name:
        cmd += f" --new-name {new_name}"
    if read_only is True:
        cmd += " --make-read-only"
    elif read_only is False:
        cmd += " --remove-read-only"

    res = await invoke(cmd)
    debug(res)

    if (not new_name) and (read_only is None):
        assert res.retcode == 1
    else:
        assert res.retcode == 0
        new_project = res.data
        if new_name:
            assert new_project["name"] == new_name
        else:
            assert new_project["name"] == name
        if read_only is True:
            assert new_project["read_only"] is True
        if read_only is False:
            assert new_project["read_only"] is False
        if read_only is None:
            assert new_project["read_only"] == project["read_only"]
