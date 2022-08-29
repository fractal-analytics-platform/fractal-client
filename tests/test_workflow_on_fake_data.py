import os
import pathlib
import shutil
import subprocess

import dask.array as da
import pytest
from devtools import debug

try:
    process = subprocess.Popen(
        ["sinfo"], stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = process.communicate()
    HAS_SLURM = True
except FileNotFoundError:
    HAS_SLURM = False


@pytest.mark.skipif(not HAS_SLURM, reason="SLURM not available")
def test_workflow_fake_data(
    # tmp_path: pathlib.Path,
):

    from fractal.fractal_cmd import dataset_update_type
    from fractal.fractal_cmd import datasets_add_resources
    from fractal.fractal_cmd import datasets_list
    from fractal.fractal_cmd import project_new
    from fractal.fractal_cmd import projects_list
    from fractal.fractal_cmd import task_add
    from fractal.fractal_cmd import task_list
    from fractal.fractal_cmd import workflow_add_task
    from fractal.fractal_cmd import workflow_apply
    from fractal.fractal_cmd import workflow_list
    from fractal.fractal_cmd import workflow_new

    testdir = os.path.dirname(__file__)
    tmp_path = pathlib.Path(f"{testdir}/tmp")

    tmp_dir = tmp_path.as_posix() + "/"
    if os.path.isdir(tmp_dir):
        shutil.rmtree(tmp_dir)

    resource_in = f"{testdir}/data/png"
    debug(tmp_dir)
    resource_out = tmp_dir

    # Quick&dirty way to ignore function decorators
    # (which are otherwise used for CLI)
    project_new = project_new.__dict__["callback"]
    projects_list = projects_list.__dict__["callback"]
    datasets_add_resources = datasets_add_resources.__dict__["callback"]
    dataset_update_type = dataset_update_type.__dict__["callback"]
    datasets_list = datasets_list.__dict__["callback"]
    task_add = task_add.__dict__["callback"]
    task_list = task_list.__dict__["callback"]
    workflow_new = workflow_new.__dict__["callback"]
    workflow_add_task = workflow_add_task.__dict__["callback"]
    workflow_list = workflow_list.__dict__["callback"]
    workflow_apply = workflow_apply.__dict__["callback"]

    # General variables and path, part 2 (do not modify this)
    project_name = "mwe-test"
    dataset_name = "dstest"
    workflow_name = "wftest"

    # Prepare and execute a workflow
    project_new(project_name, tmp_dir, dataset_name)
    projects_list()
    print()

    datasets_add_resources(project_name, dataset_name, [resource_in])
    dataset_update_type(project_name, dataset_name, "png")
    datasets_list(project_name)
    print()

    # Create workflow with only create_zarr_structure
    task_add("create_zarr_structure", "png", "zarr", "none")
    workflow_new(project_name, workflow_name, ["create_zarr_structure"])

    # Add yokogawa_to_zarr to list of tasks and to workflow
    task_add("yokogawa_to_zarr", "zarr", "zarr", "well")
    workflow_add_task(project_name, workflow_name, ["yokogawa_to_zarr"])

    workflow_list(project_name)
    print()

    workflow_apply(
        project_name,
        workflow_name,
        dataset_name,
        dataset_name,
        [resource_in],
        resource_out,
        "tests/data/parameters_workflow_on_fake_data/wf_params.json",
    )

    zarrurl = resource_out + "myplate.zarr"
    debug(zarrurl)
    assert os.path.isdir(zarrurl)

    zarrurl = resource_out + "myplate.zarr/B/03/0/0"
    data_czyx = da.from_zarr(zarrurl)
    assert data_czyx.shape == (1, 2, 2160, 2560 * 2)
    assert data_czyx[0, 0, 0, 0].compute() == 0

    shutil.rmtree(tmp_dir)


if __name__ == "__main__":
    test_workflow_fake_data()
