import os
import shutil


def test_workflow_fake_data():

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

    # General variables and paths (relative to mwe_fractal folder)
    rootdir, relative_dir = os.getcwd().split("mwe_fractal")
    if relative_dir != "":
        raise Exception(
            "ERROR: this test has hard-coded paths, it should be "
            "run from mwe_fractal folder"
        )
    resource_in = f"{rootdir}mwe_fractal/tests/data/png"
    tmp_path = f"{rootdir}mwe_fractal/tests/tmp_workflow_fake_data"
    resource_out = tmp_path

    # Remove output folder
    if os.path.isdir(tmp_path):
        print(f"Removing {tmp_path}")
        shutil.rmtree(tmp_path)

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
    project_new(project_name, tmp_path, dataset_name)
    projects_list()
    print()

    datasets_add_resources(project_name, dataset_name, [resource_in])
    dataset_update_type(project_name, dataset_name, "png")
    datasets_list(project_name)
    print()

    task_add("create_zarr_structure", "png", "zarr", "none")
    task_add("yokogawa_to_zarr", "zarr", "zarr", "well")
    workflow_new(project_name, workflow_name, ["create_zarr_structure"])
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


if __name__ == "__main__":
    test_workflow_fake_data()
