import json

import pytest


async def test_add_resource(register_user, invoke, tmp_path):

    # Create a project with its default dataset
    project_dir = str(tmp_path)
    res = await invoke(f"project new prj0 {project_dir}")
    project_id = res.data["id"]
    dataset_id = res.data["dataset_list"][0]["id"]
    assert res.data["dataset_list"][0]["resource_list"] == []

    # Add a resource
    PATH = "/some/path"
    res = await invoke(
        f"dataset add-resource {project_id} {dataset_id} {PATH}"
    )
    res.show()
    assert res.retcode == 0
    assert res.data["path"] == PATH
    assert res.data["dataset_id"] == dataset_id


async def test_add_resource_relative_path(register_user, invoke, tmp_path):
    project_dir = str(tmp_path)
    res = await invoke(f"project new prj0 {project_dir}")
    project_id = res.data["id"]
    dataset_id = res.data["dataset_list"][0]["id"]
    assert res.data["dataset_list"][0]["resource_list"] == []

    PATH = "../new/resource/path"
    with pytest.raises(ValueError):
        res = await invoke(
            f"dataset add-resource {project_id} {dataset_id} {PATH}"
        )

    PATH = "local-folder/new/resource/path"
    with pytest.raises(ValueError):
        res = await invoke(
            f"dataset add-resource {project_id} {dataset_id} {PATH}"
        )


async def test_edit_dataset(register_user, invoke, tmp_path):
    project_dir = str(tmp_path)
    res = await invoke(f"project new prj0 {project_dir}")
    project_id = res.data["id"]
    dataset_id = res.data["dataset_list"][0]["id"]

    TYPE = "this_new_type"
    NAME = "this_new_name"
    META = {"something": "else"}
    META_FILE = str(tmp_path / "meta.json")
    with open(META_FILE, "w") as f:
        json.dump(META, f)

    res = await invoke(
        f"dataset edit {project_id} {dataset_id} --new-name {NAME}"
    )
    res.show()
    assert res.data["name"] == NAME
    assert res.retcode == 0

    res = await invoke(
        f"dataset edit {project_id} {dataset_id} --new-type {TYPE}"
    )
    res.show()
    assert res.data["name"] == NAME
    assert res.data["type"] == TYPE
    assert res.retcode == 0

    res = await invoke(
        f"dataset edit {project_id} {dataset_id} --make-read-only"
    )
    res.show()
    assert res.data["name"] == NAME
    assert res.data["type"] == TYPE
    assert res.data["read_only"]
    assert res.retcode == 0

    res = await invoke(
        f"dataset edit {project_id} {dataset_id} --remove-read-only"
    )
    res.show()
    assert res.data["name"] == NAME
    assert res.data["type"] == TYPE
    assert not res.data["read_only"]
    assert res.retcode == 0

    res = await invoke(
        f"dataset edit {project_id} {dataset_id} --meta-file {META_FILE}"
    )
    res.show()
    assert res.data["name"] == NAME
    assert res.data["type"] == TYPE
    assert res.data["meta"] == META
    assert not res.data["read_only"]
    assert res.retcode == 0

    with pytest.raises(SystemExit):
        res = await invoke(
            f"dataset edit {project_id} {dataset_id} "
            "--make-read-only --remove-read-only"
        )


async def test_delete_dataset(register_user, invoke, tmp_path):
    # Create a project with its default dataset
    project_dir = str(tmp_path)
    res = await invoke(f"project new prj0 {project_dir}")
    project_id = res.data["id"]
    dataset_id = res.data["dataset_list"][0]["id"]

    # Delete dataset
    res = await invoke(f"dataset delete {project_id} {dataset_id}")

    # Check that dataset show fails
    with pytest.raises(SystemExit):
        res = await invoke(f"dataset show {project_id} {dataset_id}")


async def test_show_dataset(register_user, invoke, tmp_path):
    # Create a project with its default dataset
    project_dir = str(tmp_path)
    res = await invoke(f"project new prj0 {project_dir}")
    project_id = res.data["id"]
    dataset_id = res.data["dataset_list"][0]["id"]

    res = await invoke(f"dataset show {project_id} {dataset_id}")
    res.show()
    assert res.retcode == 0


async def test_delete_resource(register_user, invoke, tmp_path):
    # Create a project with its default dataset
    project_dir = str(tmp_path)
    res = await invoke(f"project new prj0 {project_dir}")
    project_id = res.data["id"]
    dataset_id = res.data["dataset_list"][0]["id"]
    assert res.data["dataset_list"][0]["resource_list"] == []

    # Add a resource
    PATH = "/some/path"
    res = await invoke(
        f"dataset add-resource {project_id} {dataset_id} {PATH}"
    )
    res.show()
    assert res.retcode == 0
    assert res.data["path"] == PATH
    assert res.data["dataset_id"] == dataset_id
    resource_id = res.data["id"]

    # Remove a resource
    res = await invoke(
        f"dataset rm-resource {project_id} {dataset_id} {resource_id}"
    )
    assert res.retcode == 0

    # Check that the resource was removed
    res = await invoke(f"dataset show {project_id} {dataset_id}")
    res.show()

    # Add a new resource, and check that it has the same id as the one that was
    # removed
    res = await invoke(
        f"dataset add-resource {project_id} {dataset_id} {PATH}"
    )
    res.show()
    assert res.retcode == 0
    assert resource_id == res.data["id"]
