from devtools import debug


def test_add_dataset(
    invoke,
    new_name,
    tester,
):
    # Create project
    res = invoke(f"project new {new_name()}")
    project_id = res.data["id"]

    # Create dataset (batch)
    res = invoke(
        "--batch "
        f"project add-dataset {project_id} {new_name()} "
        f"--project-dir {tester['project_dir']} --zarr-subfolder zarr"
    )
    assert res.retcode == 0
    int(res.data)

    # Create dataset (non batch)
    res = invoke(f"project add-dataset {project_id} {new_name()}")
    assert res.retcode == 0
    assert "id" in res.data.keys()


def test_edit_dataset(invoke, tester, new_name):
    res = invoke(f"project new {new_name()}")
    project_id = res.data["id"]

    res = invoke(f"project add-dataset {project_id} {new_name()}")
    dataset_id = res.data["id"]

    NAME = new_name()

    res = invoke(f"dataset edit {project_id} {dataset_id} --new-name {NAME}")
    res.show()
    assert res.data["name"] == NAME
    assert res.retcode == 0


def test_delete_dataset(invoke, new_name, tester):
    # Create a project with its default dataset
    res = invoke(f"project new {new_name()}")
    project_id = res.data["id"]

    res = invoke(f"project add-dataset {project_id} {new_name()}")
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

    res = invoke(f"project add-dataset {project_id} {new_name()}")
    dataset_id = res.data["id"]

    res = invoke(f"dataset show {project_id} {dataset_id}")
    res.show()
    assert res.retcode == 0
