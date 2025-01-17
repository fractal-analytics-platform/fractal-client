import json

from devtools import debug


def test_create_dataset(
    invoke,
    tmp_path,
    new_name,
    invoke_as_superuser,
):
    """
    Test some specific branches of the post_dataset function and parser.
    """

    res = invoke(f"project new {new_name()}")
    project_id = res.data["id"]

    attribute_filters = {"a": [1]}
    type_filters = {"b": True}

    file_attribute_filters = str(tmp_path / "attribute_filters.json")
    with open(file_attribute_filters, "w") as f:
        json.dump(attribute_filters, f)

    file_type_filters = str(tmp_path / "type_filters.json")
    with open(file_type_filters, "w") as f:
        json.dump(type_filters, f)

    res = invoke(
        (
            f"project add-dataset {project_id} {new_name()}"
            " --zarr-dir /tmp"
            f" --type-filters {file_type_filters}"
            f" --attribute-filters {file_attribute_filters}"
        )
    )

    debug(res.data)
    assert res.retcode == 0
    assert res.data["attribute_filters"] == attribute_filters
    assert res.data["type_filters"] == type_filters

    # Add a project_dir to user-settings
    res = invoke("--batch user whoami")
    assert res.retcode == 0
    user_id = res.data
    res = invoke_as_superuser(
        f"user edit {user_id} --new-project-dir /something"
    )
    assert res.retcode == 0
    res = invoke(f"--batch project add-dataset {project_id} {new_name()}")
    debug(res.data)
    assert res.retcode == 0


def test_edit_dataset(invoke, tmp_path, new_name):
    res = invoke(f"project new {new_name()}")
    project_id = res.data["id"]

    res = invoke(
        f"project add-dataset {project_id} {new_name()} --zarr-dir /tmp"
    )
    dataset_id = res.data["id"]

    NAME = new_name()

    attribute_filters = {"a": [1]}
    attribute_filters_file = str(tmp_path / "attribute_filters.json")
    with open(attribute_filters_file, "w") as f:
        json.dump(attribute_filters, f)

    type_filters = {"b": True}
    type_filters_file = str(tmp_path / "type_filters.json")
    with open(type_filters_file, "w") as f:
        json.dump(type_filters, f)

    res = invoke(f"dataset edit {project_id} {dataset_id} --new-name {NAME}")
    res.show()
    assert res.data["name"] == NAME
    assert res.retcode == 0

    res = invoke(
        f"dataset edit {project_id} {dataset_id}"
        f" --attribute-filters {attribute_filters_file}"
        f" --type-filters {type_filters_file}"
    )
    res.show()
    assert res.data["name"] == NAME
    assert res.data["type_filters"] == type_filters
    assert res.data["attribute_filters"] == attribute_filters
    assert res.retcode == 0


def test_delete_dataset(invoke, new_name):
    # Create a project with its default dataset
    res = invoke(f"project new {new_name()}")
    project_id = res.data["id"]

    res = invoke(
        f"project add-dataset {project_id} {new_name()} --zarr-dir /tmp"
    )
    dataset_id = res.data["id"]

    # Delete dataset
    res = invoke(f"dataset delete {project_id} {dataset_id}")
    debug(res.data)
    # Check that dataset show fails
    res = invoke(f"dataset show {project_id} {dataset_id}")
    assert res.data["detail"] == "Dataset not found"


def test_show_dataset(invoke, new_name):
    # Create a project with its default dataset
    res = invoke(f"project new {new_name()}")
    project_id = res.data["id"]
    res = invoke(
        f"project add-dataset {project_id} {new_name()} --zarr-dir /tmp"
    )
    dataset_id = res.data["id"]

    res = invoke(f"dataset show {project_id} {dataset_id}")
    res.show()
    assert res.retcode == 0
