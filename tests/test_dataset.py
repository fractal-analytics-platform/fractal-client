import json

from devtools import debug


def test_create_dataset(register_user, invoke, tmp_path):
    """
    Test some specific branches of the post_dataset function and parser.
    """

    res = invoke("project new prj0")
    project_id = res.data["id"]

    FILTERS = {"attributes": {"a": 1}, "types": {"b": True}}

    file_filters = str(tmp_path / "filters.json")
    with open(file_filters, "w") as f:
        json.dump(FILTERS, f)

    res = invoke(
        (
            f"project add-dataset {project_id} MyDS /tmp "
            f"--filters {file_filters}"
        )
    )
    debug(file_filters)
    debug(res.data)
    assert res.retcode == 0
    assert res.data["filters"] == FILTERS

    res = invoke(f"--batch project add-dataset {project_id} MyNewDS /tmp")
    debug(res.data)
    assert res.retcode == 0


def test_edit_dataset(register_user, invoke, tmp_path):
    res = invoke("project new prj0")
    project_id = res.data["id"]

    res = invoke(f"project add-dataset {project_id} test_name /tmp")
    dataset_id = res.data["id"]

    NAME = "this_new_name"
    FILTERS = {"attributes": {"a": 1}, "types": {"b": True}}
    FILTERS_FILE = str(tmp_path / "meta.json")
    with open(FILTERS_FILE, "w") as f:
        json.dump(FILTERS, f)

    res = invoke(f"dataset edit {project_id} {dataset_id} --new-name {NAME}")
    res.show()
    assert res.data["name"] == NAME
    assert res.retcode == 0

    res = invoke(
        f"dataset edit {project_id} {dataset_id} --filters {FILTERS_FILE}"
    )
    res.show()
    assert res.data["name"] == NAME
    assert res.data["filters"]["types"] == FILTERS["types"]
    assert res.retcode == 0


def test_delete_dataset(register_user, invoke):
    # Create a project with its default dataset
    res = invoke("project new prj0")
    project_id = res.data["id"]

    res = invoke(f"project add-dataset {project_id} test_name /tmp")
    dataset_id = res.data["id"]

    # Delete dataset
    res = invoke(f"dataset delete {project_id} {dataset_id}")
    debug(res.data)
    # Check that dataset show fails
    res = invoke(f"dataset show {project_id} {dataset_id}")
    assert res.data["detail"] == "Dataset not found"


def test_show_dataset(register_user, invoke):
    # Create a project with its default dataset
    res = invoke("project new prj0")
    project_id = res.data["id"]
    res = invoke(f"project add-dataset {project_id} test_name /tmp")
    dataset_id = res.data["id"]

    res = invoke(f"dataset show {project_id} {dataset_id}")
    res.show()
    assert res.retcode == 0
