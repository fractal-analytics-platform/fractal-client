import os

import pytest
from devtools import debug


def test_project_create(invoke, new_name):
    PROJECT_NAME = new_name()
    res = invoke(f"project new {PROJECT_NAME}")
    debug(res)
    assert res.data["name"] == PROJECT_NAME


def test_project_delete(invoke, new_name):
    # Create project
    res = invoke(f"project new {new_name()}")
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


def test_project_create_batch(invoke, new_name):
    res = invoke("project list")
    initial_projects = res.data

    res = invoke(f"--batch project new {new_name()}")
    debug(res)
    debug(res.data)
    project_id = res.data

    res = invoke("project list")
    assert len(res.data) == len(initial_projects) + 1
    assert any(project["id"] == project_id for project in res.data)


def test_project_list(invoke, new_name, tester):
    res = invoke("project list")
    initial_projects = len(res.data)

    res = invoke(f"--batch project new {new_name()}")

    zarr_dir = os.path.join(tester["project_dir"], "zarr")
    project0_id = res.data
    res = invoke(
        "--batch "
        f"project add-dataset {project0_id} {new_name()} --zarr-dir {zarr_dir}"
    )
    res = invoke(f"--batch project new {new_name()}")

    res = invoke("project list")
    debug(res)
    res.show()
    assert len(res.data) == initial_projects + 2


@pytest.mark.parametrize("patch_name", [True, False])
def test_edit_project(invoke, new_name, patch_name: bool):
    name = new_name()
    res = invoke(f"project new {name}")
    project = res.data
    project_id = project["id"]

    cmd = f"project edit {project_id}"
    if patch_name:
        NEW_NAME = new_name()
        cmd += f" --new-name {NEW_NAME}"

    res = invoke(cmd)
    debug(res)

    assert res.retcode == 0
    new_project = res.data
    if patch_name:
        assert new_project["name"] == NEW_NAME
    else:
        assert new_project["name"] == name
