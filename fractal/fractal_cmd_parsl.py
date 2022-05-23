import json
import os
from pathlib import Path
from typing import List
from typing import Optional

import click
import parsl
from devtools import debug
from parsl.addresses import address_by_hostname
from parsl.channels import LocalChannel
from parsl.config import Config
from parsl.executors import HighThroughputExecutor
from parsl.launchers import SrunLauncher
from parsl.providers import SlurmProvider
from pydantic import BaseModel
from tasks.create_zarr_structure import create_zarr_structure

# FIXME: import all from tasks (via __init__, probably)


"""
#fractal.json
{
    "fractal":{
            "version": "0.1",
            "projects":{project_name : {path: path, datasets: []}
                        },
            "tasks" : {}
            "users" : {},
            "groups": {},
             }
}

#project_name.json
{dataset_name : {resources: [],
                 type: type
                },
 "workflows" : [name : {}]
}
"""


def initialize_SlurmProvider():
    slurm = SlurmProvider(
        channel=LocalChannel(),
        nodes_per_block=1,
        cores_per_node=16,
        mem_per_node=1,  # specified in GB
        parallelism=0,
        partition="main",
        worker_init="export SQUEUE_FORMAT="
        + '"%8i %.12u %.10a %.30j %.8t %.10M %.10l %.4C %.10m %R %E"'
        + "source /opt/easybuild/software/Anaconda3/2019.07/"
        + "etc/profile.d/conda.sh\n"
        + "conda activate fractal",
        launcher=SrunLauncher(debug=True),
        walltime="01:00:00",
        cmd_timeout=60,
    )
    return slurm


def initialize_HighThroughputExecutor(provider=None, max_workers=10):
    htex = HighThroughputExecutor(
        label="htex",
        address=address_by_hostname(),
        worker_debug=True,
        max_workers=max_workers,
        provider=provider,
    )
    return htex


class TaskModel(BaseModel):
    name: str
    depends_on: Optional[str]
    input_type: str
    output_type: str


class Workflow(BaseModel):
    tasks: List[TaskModel]


def db_load():
    """
    Mocks read from global fractal database
    """
    with open("./fractal.json", "r") as f:
        db = json.load(f)
        return db


#   except FileNotFoundError:
#       return dict(projects={}, tasks={}, workflows={}, users={}, groups={})


def db_save(db):
    """
    Mocks writes to the global fractal database
    """
    with open("./fractal.json", "w") as f:
        json.dump(db, f, indent=4)


def get_project(project_name):
    """
    Loads a single project from the mock database
    """
    db = db_load()
    project_path = Path(db["fractal"]["projects"][project_name]["path"])
    ds_names = db["fractal"]["projects"][project_name]["datasets"]
    return project_path, ds_names


def project_file_load(project_name):
    """
    Loads a project file
    """
    project_path, ds_names = get_project(project_name)
    pj_file = project_name + ".json"
    with open(project_path / pj_file, "r") as f:
        prj = json.load(f)
    return prj, ds_names


def get_ds_names(project_name):
    """
    Returns datasets information from the project file
    """
    prj, ds_names = project_file_load(project_name)
    ds = prj["datasets"]
    return ds_names, ds


# def load_project_file(project_name):
#     project_path, ds_names = get_project(project_name)
#     pj_file = project_name + ".json"
#     with open(project_path / pj_file, "r") as f:
#         prj = json.load(f)
#     return prj


def save_project_file(project_name, prj_dict):
    """
    Writes project file
    """
    project_path, ds_names = get_project(project_name)
    pj_file = project_name + ".json"
    with open(project_path / pj_file, "w") as f:
        json.dump(prj_dict, f, indent=4)


def check_I_O(task_p, task_n):
    """
    Checks compatibility of tasks input/output
    """
    return task_p["output_type"] == task_n["input_type"]


@click.group()
def cli():
    pass


# ###################PROJECT#########


@cli.group()
def project():
    pass


@project.command(name="new")
@click.argument("project_name", required=True, nargs=1)
@click.argument("path", required=True, nargs=1)
@click.argument("dataset", required=True, nargs=1)
def project_new(project_name, path, dataset):
    path_obj = Path(path)
    if not path_obj.exists():
        try:
            os.mkdir(path)
        except OSError:
            print(f"Creation of the directory {path} failed")
        else:
            print(f"Successfully created the directory {path}")

    if not os.path.isfile(os.getcwd() + "/" + "fractal.json"):
        with open(os.getcwd() + "/" + "fractal.json", "w") as f:
            json.dump(
                {
                    "fractal": {
                        "version": "0.1",
                        "projects": {
                            project_name: {"path": path, "datasets": [dataset]}
                        },
                        "tasks": {},
                        "users": {},
                        "groups": {},
                    }
                },
                f,
            )

    db = db_load()

    db["fractal"]["projects"].update(
        {
            project_name: {
                "path": path_obj.resolve().as_posix(),
                "datasets": [dataset],
                "user": "",
            }
        }
    )

    pj_file = project_name + ".json"
    if not os.path.isfile(path_obj / pj_file):
        with open(path_obj / pj_file, "w") as f:
            json.dump(
                {
                    "datasets": {
                        dataset: {
                            "resources": [],
                            "type": "",
                        },
                    },
                    "workflows": {},
                },
                f,
            )

    db_save(db)


@project.command(name="list")
def projects_list():
    db = db_load()
    debug(db["fractal"]["projects"])


# #################DATASET##########


@cli.group()
def dataset():
    pass


@dataset.command(name="list")
@click.argument("project_name", required=True, nargs=1)
def datasets_list(project_name):
    ds_names, ds = get_ds_names(project_name)
    debug(ds_names, ds)


@dataset.command(name="new")
@click.argument("project_name", required=True, nargs=1)
@click.argument("dataset_name", required=True, nargs=1)
@click.argument("resources", required=True, nargs=-1)
@click.argument("ds_type", required=True, nargs=1)
def dataset_new(project_name, dataset_name, resources, ds_type):
    prj, _ = project_file_load(project_name)
    prj["datasets"].update(
        {dataset_name: {"resources": resources, "type": ds_type}}
    )
    db = db_load()
    db["fractal"]["projects"][project_name]["datasets"].append(dataset_name)
    db_save(db)
    save_project_file(project_name, prj)


@dataset.command(name="add-resources")
@click.argument("project_name", required=True, nargs=1)
@click.argument("dataset_name", required=True, nargs=1)
@click.argument("resources", required=True, nargs=-1)
def datasets_add_resources(project_name, dataset_name, resources):

    prj, _ = project_file_load(project_name)
    prj["datasets"][dataset_name]["resources"].extend(resources)
    save_project_file(project_name, prj)


@dataset.command(name="update-type")
@click.argument("project_name", required=True, nargs=1)
@click.argument("dataset_name", required=True, nargs=1)
@click.argument("ds_type", required=True, nargs=1)
def dataset_update_type(project_name, dataset_name, ds_type):

    prj, _ = project_file_load(project_name)
    prj["datasets"][dataset_name]["type"] = ds_type
    save_project_file(project_name, prj)


# ################TASK##################
@cli.group()
def task():
    pass


@task.command(name="list")
def task_list():
    db = db_load()
    debug(db["fractal"]["tasks"])


@task.command(name="add")
@click.argument("task_name", required=True, nargs=1)
@click.argument("input_type", required=True, nargs=1)
@click.argument("output_type", required=True, nargs=1)
@click.argument("depends_on", required=False)
def task_add(task_name, input_type, output_type, depends_on=None):
    db = db_load()
    db["fractal"]["tasks"][task_name] = TaskModel(
        name=task_name,
        input_type=input_type,
        output_type=output_type,
        depends_on=depends_on,
        # FIXME: path/function...
        # alll we need to obtain the *python function* to call
    ).dict()
    # FIXME WRITE IN tasks/dict_tasks.py
    # FIXME: re-import tasks/dict_tasks.py
    db_save(db)


# ########################WORKFLOW######


@cli.group()
def workflow():
    pass


@workflow.command(name="list")
@click.argument("project_name", required=True, nargs=1)
def workflow_list(project_name):
    prj, _ = project_file_load(project_name)
    debug(prj["workflows"])


@workflow.command(name="new")
@click.argument("project_name", required=True, nargs=1)
@click.argument("workflow_name", required=True, nargs=1)
@click.argument("tasks", required=True, nargs=-1)
def workflow_new(project_name, workflow_name, tasks):

    db = db_load()
    if len(tasks) > 1:
        task_p = db["fractal"]["tasks"][tasks[0]]
        for task_n in db["fractal"]["tasks"][tasks[1:]]:
            if not check_I_O(task_p, task_n):
                raise
            task_p = task_n

    prj, _ = project_file_load(project_name)

    prj["workflows"].update(
        {
            workflow_name: Workflow(
                tasks=[
                    db["fractal"]["tasks"][task_name] for task_name in tasks
                ],
            ).dict()
        }
    )
    save_project_file(project_name, prj)


@workflow.command(name="add-task")
@click.argument("project_name", required=True, nargs=1)
@click.argument("workflow_name", required=True, nargs=1)
@click.argument("tasks", required=True, nargs=-1)
def workflow_add_task(project_name, workflow_name, tasks):

    db = db_load()
    prj, _ = project_file_load(project_name)

    if len(tasks) > 1:
        task_p = db["fractal"]["tasks"][tasks[0]]
        for task in tasks[1:]:
            task_n = db["fractal"]["tasks"][task]
            if not check_I_O(task_p, task_n):
                raise
            task_p = task_n
            prj["workflows"][workflow_name]["tasks"].append(task_p)

    task = db["fractal"]["tasks"][tasks[0]]
    prj["workflows"][workflow_name]["tasks"].append(task)

    save_project_file(project_name, prj)


@workflow.command(name="apply")
@click.argument("project_name", required=True, nargs=1)
@click.argument("workflow_name", required=True, nargs=1)
@click.argument("input_dataset", required=True, nargs=1)
@click.argument("output_dataset", required=True, nargs=1)
@click.argument("resource_in", required=True, nargs=1)
@click.argument("resource_out", required=True, nargs=1)
def workflow_apply(
    project_name,
    workflow_name,
    input_dataset,
    output_dataset,
    resource_in,
    resource_out,
):

    prj, _ = project_file_load(project_name)

    # Verify that resource_in has been added to the resources of input_dataset
    dataset_resources = prj["datasets"][input_dataset]["resources"]
    if resource_in not in dataset_resources:
        raise

    # If the resource_out folder is not there, create it
    path_resource_out = Path(resource_out)
    if not path_resource_out.exists():
        os.mkdir(path_resource_out)
        prj["datasets"][output_dataset]["resources"].extend(
            [path_resource_out.resolve().as_posix()]
        )
        save_project_file(project_name, prj)

    # Collect tasks
    task_names = [t["name"] for t in prj["workflows"][workflow_name]["tasks"]]
    debug(task_names)

    # Hard-coded parameters
    ext = prj["datasets"][input_dataset]["type"]
    num_levels = 5
    dims = (1, 1)
    coarsening_factor_xy = 2
    coarsening_factor_z = 1

    # Task 0: create zarr
    zarrurl_list, chl_list = create_zarr_structure(
        in_path=resource_in,
        out_path=resource_out,
        ext=ext,
        num_levels=num_levels,
        dims=dims,
    )
    print(f"zarrurl_list: {zarrurl_list}")

    # FIXME validate tasks somewhere

    # FIXME: define slurm provider
    provider = initialize_SlurmProvider()
    htex = initialize_HighThroughputExecutor(provider=provider)
    config = Config(executors=[htex])
    parsl.clear()
    parsl.load(config)

    # Loop over tasks
    # db = db_load()
    for task in task_names:
        futures = []

        if task == "yokogawa_to_zarr":
            # FIXME: add task-specific wrapper
            # QUESTION: where will parsl act? probably within the wrapper??
            kwargs = dict(
                in_path=resource_in,
                out_path=resource_out,
                ext=ext,
                dims=dims,
                chl_list=chl_list,
                coarsening_factor_xy=coarsening_factor_xy,
                coarsening_factor_z=coarsening_factor_z,
            )
            # function = yokogawa_to_zarr   # <-- how to do it more generally??

            def function(zarrurl, **kwargs_):
                from fractal.tasks.yokogawa_to_zarr import yokogawa_to_zarr

                return yokogawa_to_zarr(zarrurl, **kwargs_)

            # function = DICT_TASKS[task]  #FIXME
        # elif task == "illumination":
        #    kwargs = dict(..
        #                  threshold=0.5)
        #    function = illumination

        # app1 = parsl.python_app(function=task_increment)(1)
        # app2 = PythonApp(task_increment)(2)

        for zarrurl in zarrurl_list:
            future = parsl.python_app(function=function)(zarrurl, **kwargs)
            # future = function(zarrurl, **kwargs)
            futures.append(future)

        print(futures)
        [app.result() for app in futures]
        print(futures)
        print()


@cli.group()
def ps():
    pass


@ps.command(name="list")
@click.argument("workflow_name", required=True, nargs=1)
def ps_list(workflow_name):
    raise NotImplementedError


if __name__ == "__main__":
    cli()
