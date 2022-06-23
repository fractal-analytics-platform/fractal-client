import json
import os
from pathlib import Path
from typing import List
from typing import Optional

import click
import parsl
from devtools import debug
from parsl import python_app
from parsl.config import Config
from parsl_config import define_HighThroughputExecutor
from parsl_config import define_MonitoringHub
from parsl_config import define_SlurmProvider
from pydantic import BaseModel

import fractal.fractal_config as fractal_config


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


def add_slash_to_path(path):
    if path.endswith("/"):
        return path
    else:
        return path + "/"


class TaskModel(BaseModel):
    name: str
    depends_on: Optional[str]
    input_type: str
    output_type: str
    parallelization_level: str


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

    path = add_slash_to_path(path)

    path_obj = Path(path)
    if not path_obj.exists():
        try:
            os.mkdir(path)
        except OSError:
            print(f"Creation of the directory {path} failed")
        else:
            print(f"Successfully created the directory {path}")

    if not os.path.isfile(
        os.getcwd() + "/" + "fractal.json"
    ):  # FIXME move it to path??
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

    with open("dictionary_tasks.py", "w") as out:
        out.write("dict_tasks = {}\n\n")


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

    sanitized_resources = [add_slash_to_path(res) for res in resources]

    prj, _ = project_file_load(project_name)
    prj["datasets"].update(
        {dataset_name: {"resources": sanitized_resources, "type": ds_type}}
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

    sanitized_resources = [add_slash_to_path(res) for res in resources]

    prj, _ = project_file_load(project_name)
    prj["datasets"][dataset_name]["resources"].extend(sanitized_resources)
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
@click.argument("parallelization_level", required=True, nargs=1)
@click.argument("depends_on", required=False)
def task_add(
    task_name, input_type, output_type, parallelization_level, depends_on=None
):
    db = db_load()
    db["fractal"]["tasks"][task_name] = TaskModel(
        name=task_name,
        input_type=input_type,
        output_type=output_type,
        depends_on=depends_on,
        parallelization_level=parallelization_level,
    ).dict()

    with open("dictionary_tasks.py", "a") as out:  # FIXME add path
        out.write(f"from fractal.tasks.{task_name} import {task_name}\n")
        out.write(f'dict_tasks["{task_name}"] = {task_name}\n')
        out.write("\n")

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
                raise Exception("I/O error in workflow_new")
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
                raise Exception("I/O error in workflow_add_task")
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
@click.argument("resources_in", required=True, nargs=-1)
@click.argument("resource_out", required=True, nargs=1)
@click.argument("json_worker_params", required=True, nargs=1)
def workflow_apply(
    project_name,
    workflow_name,
    input_dataset,
    output_dataset,
    resources_in,
    resource_out,
    json_worker_params,
):

    resources_in = [
        add_slash_to_path(resource_in) for resource_in in resources_in
    ]
    resource_out = add_slash_to_path(resource_out)

    prj, _ = project_file_load(project_name)

    # Verify that resources_in has been added to the resources of input_dataset
    dataset_resources = prj["datasets"][input_dataset]["resources"]
    for resource_in in resources_in:
        if resource_in not in dataset_resources:
            raise Exception(
                f"Error in workflow_apply, {resource_in} not in"
                f" {dataset_resources}"
            )

    # If the resource_out folder is not there, create it
    path_resource_out = Path(resource_out)
    if not path_resource_out.exists():
        os.mkdir(path_resource_out)
        prj["datasets"][output_dataset]["resources"].extend(
            [add_slash_to_path(path_resource_out.resolve().as_posix())]
        )
        save_project_file(project_name, prj)

    # Collect tasks
    task_names = [t["name"] for t in prj["workflows"][workflow_name]["tasks"]]
    debug(task_names)

    # Hard-coded parameters
    with open(json_worker_params, "r") as file_params:
        params = json.load(file_params)
    ext = prj["datasets"][input_dataset]["type"]
    coarsening_xy = params["coarsening_xy"]
    coarsening_z = params["coarsening_z"]
    num_levels = params["num_levels"]
    rows, cols = params["dims"]
    workflow_name = params["workflow_name"]
    path_dict_channels = params["channel_file"]
    path_dict_corr = params["path_dict_corr"]
    delete_input = params.get("delete_input", False)

    # FIXME validate tasks somewhere?

    # Hard-coded parsl options
    OPENBLAS_NUM_THREADS = "1"  # FIXME
    fmt = "%8i %.12u %.10a %.30j %.8t %.10M %.10l %.4C %.10m %R %E"
    os.environ["SQUEUE_FORMAT"] = fmt

    # fractal_config has defined worker_init and some slurm options
    provider = define_SlurmProvider(
        nodes_per_block=fractal_config.nodes_per_block,
        cores_per_node=fractal_config.cores_per_node,
        mem_per_node_GB=fractal_config.mem_per_node_GB,
        partition=fractal_config.partition,
        worker_init=fractal_config.worker_init,
        max_blocks=fractal_config.max_blocks,
        exclusive=fractal_config.exclusive,
    )
    htex = define_HighThroughputExecutor(
        provider=provider, max_workers=fractal_config.max_workers
    )
    monitoring = define_MonitoringHub(workflow_name=workflow_name)
    config = Config(executors=[htex], monitoring=monitoring)
    # config = Config(executors=[htex])
    parsl.clear()
    parsl.load(config)

    from fractal.dictionary_tasks import dict_tasks

    debug(dict_tasks)

    @parsl.python_app
    def collect_intermediate_results(inputs=[]):
        return [x for x in inputs]

    # Task 0
    if task_names[0] in [
        "create_zarr_structure",
        "create_zarr_structure_multifov",
    ]:
        kwargs = dict(
            in_paths=resources_in,
            out_path=resource_out,
            path_dict_channels=path_dict_channels,
            ext=ext,
            num_levels=num_levels,
        )

        @parsl.python_app
        def app_create_zarr_structure(**kwargs_):
            os.environ["OPENBLAS_NUM_THREADS"] = OPENBLAS_NUM_THREADS
            import fractal.dictionary_tasks  # noqa: F401

            return dict_tasks[task_names[0]](**kwargs)

        future = app_create_zarr_structure(**kwargs)

        if task_names[0] == "create_zarr_structure":
            zarrurls, chl_list = future.result()
            debug(zarrurls)
            debug(chl_list)
        elif task_names[0] == "create_zarr_structure_multifov":
            zarrurls, chl_list, well_to_sites = future.result()
            debug(zarrurls)
            debug(chl_list)
            debug(well_to_sites)
        task_names = task_names[1:]  # FIXME
    else:
        print(
            "ERROR/WARNING: Workflows must start with create_zarr_structure*"
        )

    # Tasks 1,2,...
    db = db_load()
    for task in task_names:
        if len(resources_in) > 1:
            raise Exception(
                "ERROR\n"
                "Support for len(resources_in)>1 is not there.\n"
                "Hint: we should modify the in_path argument of"
                "yokogawa_to_zarr and take a global dict of paths."
            )

        if task == "yokogawa_to_zarr":
            kwargs = dict(
                in_path=resources_in[0],  # FIXME
                ext=ext,
                delete_input=delete_input,
                rows=rows,
                cols=cols,
                chl_list=chl_list,
                num_levels=num_levels,
                coarsening_xy=coarsening_xy,
                coarsening_z=coarsening_z,
            )
        if task == "yokogawa_to_zarr_multifov":
            kwargs = dict(
                in_path=resources_in[0],  # FIXME
                ext=ext,
                delete_input=delete_input,
                chl_list=chl_list,
                # sites_dict=well_to_sites,
                num_levels=num_levels,
                coarsening_xy=coarsening_xy,
                coarsening_z=coarsening_z,
            )

        elif task == "maximum_intensity_projection":
            kwargs = dict(
                coarsening_xy=coarsening_xy,
            )
        elif task == "replicate_zarr_structure_mip":
            kwargs = {}
        elif task == "replicate_zarr_structure":
            kwargs = dict(newzarrurl="new")
        elif task == "illumination_correction":
            kwargs = dict(
                chl_list=chl_list,
                path_dict_corr=path_dict_corr,
                coarsening_xy=coarsening_xy,
                overwrite=True,
                # background=background,
            )

        @python_app
        def app(zarrurl, **kwargs_):
            os.environ["OPENBLAS_NUM_THREADS"] = OPENBLAS_NUM_THREADS
            import fractal.dictionary_tasks  # noqa: F401

            return dict_tasks[task](zarrurl, **kwargs_)

        futures = []
        parallelization_level = db["fractal"]["tasks"][task][
            "parallelization_level"
        ]
        for zarrurl in zarrurls[parallelization_level]:
            future = app(zarrurl, **kwargs)
            futures.append(future)

        print(futures)
        # [future.result() for future in futures]
        collect_intermediate_results(inputs=futures).result()
        print(futures)
        print()

        # FIXME task-specific "naming" of output
        # FIXME to be validated

    # FIXME: Is this needed??
    parsl.clear()


@cli.group()
def ps():
    pass


@ps.command(name="list")
@click.argument("workflow_name", required=True, nargs=1)
def ps_list(workflow_name):
    raise NotImplementedError


if __name__ == "__main__":
    cli()
