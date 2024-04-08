import pytest
from devtools import debug


def test_project_create(register_user, invoke):
    PROJECT_NAME = "project_name"
    res = invoke(f"project new {PROJECT_NAME}")
    debug(res)
    assert res.data["name"] == PROJECT_NAME


def test_project_delete(register_user, invoke):

    # Create project
    res = invoke("project new MyProj1")
    res.show()
    project_id_1 = res.data["id"]

    # Show project
    res = invoke(f"project show {project_id_1}")
    res.show()

    # Delete project
    res = invoke(f"project delete {project_id_1}")
    assert res.retcode == 0

    # Try to show deleted project, and fail
    with pytest.raises(SystemExit):
        res = invoke(f"project show {project_id_1}")


def test_project_create_batch(register_user, invoke):
    res = invoke("--batch project new MyProj1")
    debug(res)
    debug(res.data)
    project_id = int(res.data)
    assert project_id == 1


def test_project_list(register_user, invoke):
    res = invoke("project list")
    debug(res)
    assert len(res.data) == 0

    res.show()

    res = invoke("--batch project new proj0")
    project0_id = res.data
    res = invoke(f"--batch project add-dataset {project0_id} NAME /tmp")
    res = invoke("--batch project new proj1")

    res = invoke("project list")
    debug(res)
    res.show()
    assert len(res.data) == 2


@pytest.mark.parametrize("new_name", ["new_name", None])
def test_edit_project(
    register_user,
    invoke,
    new_name,
    tmp_path,
):
    name = "name"
    res = invoke(f"project new {name}")
    project = res.data
    project_id = project["id"]

    cmd = f"project edit {project_id}"
    if new_name:
        cmd += f" --new-name {new_name}"

    res = invoke(cmd)
    debug(res)

    assert res.retcode == 0
    new_project = res.data
    if new_name:
        assert new_project["name"] == new_name
    else:
        assert new_project["name"] == name
