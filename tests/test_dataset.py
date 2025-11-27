import os

from devtools import debug


def test_create_dataset(
    invoke,
    new_name,
    invoke_as_superuser,
    tester,
):
    """
    Test some specific branches of the post_dataset function and parser.
    """

    res = invoke(f"project new {new_name()}")
    project_id = res.data["id"]

    zarr_dir = os.path.join(tester["project_dir"], "zarr")
    res = invoke(
        f"project add-dataset {project_id} {new_name()} --zarr-dir {zarr_dir}"
    )

    debug(res.data)
    assert res.retcode == 0

    # Add a project_dir to user
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


def test_edit_dataset(invoke, tester, new_name):
    res = invoke(f"project new {new_name()}")
    project_id = res.data["id"]

    zarr_dir = os.path.join(tester["project_dir"], "zarr")
    res = invoke(
        f"project add-dataset {project_id} {new_name()} --zarr-dir {zarr_dir}"
    )
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

    zarr_dir = os.path.join(tester["project_dir"], "zarr")

    res = invoke(
        f"project add-dataset {project_id} {new_name()} --zarr-dir {zarr_dir}"
    )
    dataset_id = res.data["id"]

    # Delete dataset
    res = invoke(f"dataset delete {project_id} {dataset_id}")
    debug(res.data)
    # Check that dataset show fails
    res = invoke(f"dataset show {project_id} {dataset_id}")
    assert res.data["detail"] == "Dataset not found"


def test_show_dataset(invoke, new_name, tester):
    # Create a project with its default dataset
    res = invoke(f"project new {new_name()}")
    project_id = res.data["id"]

    zarr_dir = os.path.join(tester["project_dir"], "zarr")
    res = invoke(
        f"project add-dataset {project_id} {new_name()} --zarr-dir {zarr_dir}"
    )
    dataset_id = res.data["id"]

    res = invoke(f"dataset show {project_id} {dataset_id}")
    res.show()
    assert res.retcode == 0
