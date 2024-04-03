import json

import pytest
from devtools import debug


def test_create_dataset(register_user, invoke, tmp_path):
    """
    Test some specific branches of the post_dataset function and parser.
    """

    res = invoke("project new prj0")
    project_id = res.data["id"]

    METADATA = {"some": "value"}
    TYPE = "sometype"

    file_metadata = str(tmp_path / "metadata.json")
    with open(file_metadata, "w") as f:
        json.dump(METADATA, f)

    res = invoke(
        (
            f"project add-dataset {project_id} MyDS "
            f"--metadata {file_metadata} "
            f"--type {TYPE}"
        )
    )
    debug(res.data)
    assert res.retcode == 0
    assert res.data["meta"] == METADATA
    assert res.data["type"] == TYPE

    res = invoke(f"--batch project add-dataset {project_id} MyNewDS")
    debug(res.data)
    assert res.retcode == 0


def test_add_resource(register_user, invoke):

    # Create a project with its default dataset
    res = invoke("project new prj0")
    project_id = res.data["id"]

    res = invoke(f"project add-dataset {project_id} test_name")
    dataset_id = res.data["id"]
    assert res.data["resource_list"] == []

    # Add a resource
    PATH = "/some/path"
    res = invoke(f"dataset add-resource {project_id} {dataset_id} {PATH}")
    res.show()
    assert res.retcode == 0
    assert res.data["path"] == PATH
    assert res.data["dataset_id"] == dataset_id


def test_add_resource_relative_path(register_user, invoke):
    res = invoke("project new prj0")
    project_id = res.data["id"]

    res = invoke(f"project add-dataset {project_id} test_name")
    dataset_id = res.data["id"]
    assert res.data["resource_list"] == []

    PATH = "../new/resource/path"
    with pytest.raises(SystemExit):
        res = invoke(f"dataset add-resource {project_id} {dataset_id} {PATH}")

    PATH = "local-folder/new/resource/path"
    with pytest.raises(SystemExit):
        res = invoke(f"dataset add-resource {project_id} {dataset_id} {PATH}")


def test_edit_dataset(register_user, invoke, tmp_path):
    res = invoke("project new prj0")
    project_id = res.data["id"]

    res = invoke(f"project add-dataset {project_id} test_name")
    dataset_id = res.data["id"]

    TYPE = "this_new_type"
    NAME = "this_new_name"
    META = {"something": "else"}
    META_FILE = str(tmp_path / "meta.json")
    with open(META_FILE, "w") as f:
        json.dump(META, f)

    res = invoke(f"dataset edit {project_id} {dataset_id} --new-name {NAME}")
    res.show()
    assert res.data["name"] == NAME
    assert res.retcode == 0

    res = invoke(f"dataset edit {project_id} {dataset_id} --new-type {TYPE}")
    res.show()
    assert res.data["name"] == NAME
    assert res.data["type"] == TYPE
    assert res.retcode == 0

    res = invoke(f"dataset edit {project_id} {dataset_id} --make-read-only")
    res.show()
    assert res.data["name"] == NAME
    assert res.data["type"] == TYPE
    assert res.data["read_only"]
    assert res.retcode == 0

    res = invoke(f"dataset edit {project_id} {dataset_id} --remove-read-only")
    res.show()
    assert res.data["name"] == NAME
    assert res.data["type"] == TYPE
    assert not res.data["read_only"]
    assert res.retcode == 0

    res = invoke(
        f"dataset edit {project_id} {dataset_id} --meta-file {META_FILE}"
    )
    res.show()
    assert res.data["name"] == NAME
    assert res.data["type"] == TYPE
    assert res.data["meta"] == META
    assert not res.data["read_only"]
    assert res.retcode == 0

    with pytest.raises(SystemExit):
        res = invoke(
            f"dataset edit {project_id} {dataset_id} "
            "--make-read-only --remove-read-only"
        )


def test_delete_dataset(register_user, invoke):
    # Create a project with its default dataset
    res = invoke("project new prj0")
    project_id = res.data["id"]

    res = invoke(f"project add-dataset {project_id} test_name")
    dataset_id = res.data["id"]

    # Delete dataset
    res = invoke(f"dataset delete {project_id} {dataset_id}")
    debug(res.data)
    # Check that dataset show fails
    # with pytest.raises(SystemExit):
    res = invoke(f"dataset show {project_id} {dataset_id}")
    assert res.data["detail"] == 'Dataset not found'


def test_show_dataset(register_user, invoke):
    # Create a project with its default dataset
    res = invoke("project new prj0")
    project_id = res.data["id"]
    res = invoke(f"project add-dataset {project_id} test_name")
    dataset_id = res.data["id"]

    res = invoke(f"dataset show {project_id} {dataset_id}")
    res.show()
    assert res.retcode == 0


def test_delete_resource(register_user, invoke):
    # Create a project with its default dataset
    res = invoke("project new prj0")
    project_id = res.data["id"]
    res = invoke(f"project add-dataset {project_id} test_name")
    dataset_id = res.data["id"]
    assert res.data["resource_list"] == []

    # Add a resource
    PATH = "/some/path"
    res = invoke(f"dataset add-resource {project_id} {dataset_id} {PATH}")
    res.show()
    assert res.retcode == 0
    assert res.data["path"] == PATH
    assert res.data["dataset_id"] == dataset_id
    resource_id = res.data["id"]

    # Show dataset
    res = invoke(f"dataset show {project_id} {dataset_id}")
    res.show()

    # Remove a resource
    res = invoke(
        f"dataset rm-resource {project_id} {dataset_id} {resource_id}"
    )
    assert res.retcode == 0

    # Check that the resource was removed
    res = invoke(f"dataset show {project_id} {dataset_id}")
    res.show()

    # Add a new resource, and check that it has the same id as the one that was
    # removed
    res = invoke(
        f"--batch dataset add-resource {project_id} {dataset_id} {PATH}"
    )
    assert res.retcode == 0
    assert res.data == resource_id


def test_dataset_history_command(register_user, invoke):
    """
    Only test the client interface, not the fractal-server business logic.
    """
    res = invoke("project new prj0")
    project_id = res.data["id"]
    res = invoke(f"project add-dataset {project_id} test_name")
    dataset_id = res.data["id"]

    debug(f"dataset history {project_id} {dataset_id}")
    res = invoke(f"dataset history {project_id} {dataset_id}")
    res.show()
    assert res.retcode == 0


def test_dataset_status_command(register_user, invoke):
    """
    Only test the client interface, not the fractal-server business logic.
    """
    res = invoke("project new prj0")
    project_id = res.data["id"]
    res = invoke(f"project add-dataset {project_id} test_name")
    dataset_id = res.data["id"]

    res = invoke(f"dataset status {project_id} {dataset_id}")
    res.show()
    assert res.retcode == 0
