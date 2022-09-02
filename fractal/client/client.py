"""
Copyright 2022 (C) Friedrich Miescher Institute for Biomedical Research and
University of Zurich

Original authors:
Jacopo Nespolo <jacopo.nespolo@exact-lab.it>
Marco Franzon <marco.franzon@exact-lab.it>
Tommaso Comparin <tommaso.comparin@exact-lab.it>

This file is part of Fractal and was originally developed by eXact lab S.r.l.
<exact-lab.it> under contract with Liberali Lab from the Friedrich Miescher
Institute for Biomedical Research and Pelkmans Lab from the University of
Zurich.
"""
import json
from typing import Any
from typing import Dict
from typing import List
from typing import Optional

import asyncclick as click
import httpx
from rich import print_json
from rich.console import Console
from rich.table import Table

from ._auth import AuthToken
from .config import settings
from fractal.common.models import ResourceRead
from fractal.common.models import SubtaskCreate
import logging

console = Console()


async def _extract_project_and_dataset(project_name: str, dataset_name: str):
    async with httpx.AsyncClient() as client:
        auth = AuthToken(client=client)
        res = await client.get(
            f"{settings.BASE_URL}/project/",
            headers=await auth.header(),
        )
        projects = res.json()
        # Find project
        try:
            project = [p for p in projects if p["name"] == project_name][0]
            project_id = project["id"]
        except IndexError as e:
            raise IndexError(f"Project {project_name} not found", str(e))
        # Get dataset list
        dataset_list = [
            project["dataset_list"]
            for project in projects
            if project["id"] == project_id
        ][0]
        # Find dataset
        try:
            dataset = [
                ds for ds in dataset_list if ds["name"] == dataset_name
            ][0]
        except IndexError as e:
            raise IndexError(f"Dataset {dataset_name} not found", str(e))

    return project, dataset


@click.group()
async def cli():
    pass


@cli.command(name="login")
async def login():
    async with httpx.AsyncClient() as client:
        auth = AuthToken(client=client)
        await auth()
        logging.debug(await auth.header())


# PROJECT GROUP


@cli.group()
async def project():
    pass


@project.command(name="new")
@click.argument("name", required=True, nargs=1)
@click.argument(
    "path",
    required=True,
    nargs=1,
)
@click.option(
    "--dataset",
    required=True,
    nargs=1,
    default="default",
    help=(
        "name of first dataset. By default, the dataset `default` is "
        "created and added to the project"
    ),
)
async def project_new(name: str, path: str, dataset: str) -> None:
    """
    Create new project, together with its first dataset

    NAME (str): project name

    PATH (str): project path, i.e., the path where all the artifacts will be
    saved.
    """
    from fractal.common.models import ProjectBase

    project = ProjectBase(name=name, project_dir=path)
    async with httpx.AsyncClient() as client:
        auth = AuthToken(client=client)
        res = await client.post(
            f"{settings.BASE_URL}/project/",
            json=project.dict(),
            headers=await auth.header(),
        )
        logging.debug(res.status_code)
        if res.status_code != 201:
            raise Exception(
                "ERROR (hint: maybe the project already exists?)", res
            )
        print_json(data=res.json())


@project.command(name="list")
async def project_list():
    from fractal.common.models import ProjectRead

    async with httpx.AsyncClient() as client:
        auth = AuthToken(client=client)
        res = await client.get(
            f"{settings.BASE_URL}/project/",
            headers=await auth.header(),
        )
        data = res.json()
        project_list = [ProjectRead(**item) for item in data]

        table = Table(title="Project List")
        table.add_column("Id", style="cyan", no_wrap=True)
        table.add_column("Name", style="magenta")
        table.add_column("Proj. Dir.", justify="right", style="green")
        table.add_column("Dataset list", style="white")
        table.add_column("Read only", justify="center")

        for p in project_list:
            # Map p.read_only (True/False) to read_only_icon (✅/❌)
            if p.read_only:
                read_only_icon = "✅"
            else:
                read_only_icon = "❌"

            p_dataset_list = str([dataset.name for dataset in p.dataset_list])

            table.add_row(
                str(p.id),
                p.name,
                p.project_dir,
                str(p_dataset_list),
                read_only_icon,
            )

        console.print(table)


@project.command(name="add-dataset")
@click.argument("project_name", required=True, nargs=1, type=str)
@click.argument(
    "dataset_name",
    type=str,
    required=True,
    nargs=1,
)
@click.option(
    "--meta",
    nargs=1,
    default=None,
    help="JSON file with meta",
)
@click.option(
    "--type",
    nargs=1,
    default="zarr",
    help=("The type of objects into the dataset"),
)
async def add_dataset(
    project_name: str,
    dataset_name: str,
    meta: Dict[str, Any],
    type: Optional[str],
) -> None:
    """
    Add an existing dataset to an existing project
    """
    from fractal.common.models import DatasetCreate

    if meta is None:
        meta_json = {}
    else:
        with open(meta, "r", encoding="utf-8") as json_file:
            meta_json = json.load(json_file)

    # Find project_id
    async with httpx.AsyncClient() as client:
        auth = AuthToken(client=client)
        res = await client.get(
            f"{settings.BASE_URL}/project/",
            headers=await auth.header(),
        )
        projects = res.json()
        # Find project
        try:
            project = [p for p in projects if p["name"] == project_name][0]
            project_id = project["id"]
        except IndexError as e:
            raise IndexError(f"Project {project_name} not found", str(e))

    # Check that there is no other dataset with the same name
    from fractal.common.models import ProjectRead

    async with httpx.AsyncClient() as client:
        auth = AuthToken(client=client)
        res = await client.get(
            f"{settings.BASE_URL}/project/",
            headers=await auth.header(),
        )
        data = res.json()
        project_list = [ProjectRead(**item) for item in data]
        logging.debug(project_list)
        logging.debug(project_id)
        try:
            project = [p for p in project_list if p.id == project_id][0]
        except IndexError as e:
            raise IndexError("No project found with this project_id", str(e))
        list_dataset_names = [dataset.name for dataset in project.dataset_list]

        if dataset_name in list_dataset_names:
            raise Exception(
                f"Dataset name {dataset_name} already in use, "
                "pick another one"
            )

    dataset = DatasetCreate(
        name=dataset_name, project_id=project_id, type=type, meta=meta_json
    )
    async with httpx.AsyncClient() as client:
        auth = AuthToken(client=client)
        res = await client.post(
            f"{settings.BASE_URL}/project/{project_id}/",
            json=dataset.dict(),
            headers=await auth.header(),
        )
        print_json(data=res.json())


# DATASET GROUP


@cli.group()
async def dataset():
    pass


@dataset.command(name="show")
@click.argument("project_name", required=True, type=str, nargs=1)
@click.argument("dataset_name", required=True, type=str, nargs=1)
async def dataset_show(project_name: str, dataset_name: str) -> None:
    """
    Show details of an existing dataset
    """

    project, dataset = await _extract_project_and_dataset(
        project_name, dataset_name
    )

    table = Table(title="Dataset")
    table.add_column("Id", style="cyan", no_wrap=True)
    table.add_column("Name", justify="right", style="green")
    table.add_column("Type", style="white")
    table.add_column("Meta", justify="center")
    table.add_column("Read only", justify="center")

    if dataset["read_only"]:
        ds_read_only = "✅"
    else:
        ds_read_only = "❌"

    table.add_row(
        str(dataset["id"]),
        dataset["name"],
        dataset["type"],
        str(dataset["meta"]),
        ds_read_only,
    )

    console.print(table)


@dataset.command(name="add-resource")
@click.argument("project_name", required=True, type=str, nargs=1)
@click.argument("dataset_name", required=True, type=str, nargs=1)
@click.argument("path", required=True, type=str, nargs=1)
@click.option("--glob_pattern", required=True, type=str, default="", nargs=1)
async def add_resource(
    project_name: str, dataset_name: str, path: str, glob_pattern: str
):

    """
    Add a new resource to an exisisting dataset
    """

    from fractal.common.models import ResourceCreate

    resource = ResourceCreate(path=path, glob_pattern=glob_pattern)

    project, dataset = await _extract_project_and_dataset(
        project_name, dataset_name
    )
    project_id = project["id"]
    dataset_id = dataset["id"]

    async with httpx.AsyncClient() as client:
        auth = AuthToken(client=client)
        res = await client.post(
            f"{settings.BASE_URL}/project/{project_id}/{dataset_id}",
            json=resource.dict(),
            headers=await auth.header(),
        )

        print_json(data=res.json())


@dataset.command(name="show-resources")
@click.argument("project_name", required=True, type=str, nargs=1)
@click.argument("dataset_name", required=True, type=str, nargs=1)
async def get_resource(
    project_name: str,
    dataset_name: str,
):

    project, dataset = await _extract_project_and_dataset(
        project_name, dataset_name
    )
    project_id = project["id"]
    dataset_id = dataset["id"]

    async with httpx.AsyncClient() as client:
        auth = AuthToken(client=client)
        res = await client.get(
            f"{settings.BASE_URL}/project/{project_id}/{dataset_id}",
            headers=await auth.header(),
        )
        data = res.json()
        resource_list = [ResourceRead(**item) for item in data]

        table = Table(title="Resource List")
        table.add_column("Id", style="cyan", no_wrap=True)
        table.add_column("Path", justify="right", style="green")
        table.add_column("Dataset Id", style="white")
        table.add_column("Glob Pattern", style="red")
        for r in resource_list:

            table.add_row(str(r.id), r.path, str(r.dataset_id), r.glob_pattern)

        console.print(table)


@dataset.command(name="modify-dataset")
@click.argument("project_name", required=True, type=str, nargs=1)
@click.argument("dataset_name", required=True, type=str, nargs=1)
@click.option(
    "--new_dataset_name",
    nargs=1,
)
@click.option(
    "--meta",
    nargs=1,
)
@click.option(
    "--type",
    nargs=1,
    help=("The type of objects into the dataset"),
)
@click.option(
    "--read_only",
    nargs=1,
    help=("Writing permissions"),
)
async def modify_dataset(
    project_name: str,
    dataset_name: str,
    new_dataset_name: str = None,
    meta: str = None,
    type: str = None,
    read_only: bool = None,
):

    if meta is None:
        mt = None
    else:
        with open(meta, "r", encoding="utf-8") as m:
            mt = json.load(m)

    updates = dict(
        name=new_dataset_name,
        meta=mt,
        type=type,
        read_only=read_only,
    )
    updates_not_none = {
        key: value for key, value in updates.items() if value is not None
    }

    project, dataset = await _extract_project_and_dataset(
        project_name, dataset_name
    )
    project_id = project["id"]
    dataset_id = dataset["id"]

    async with httpx.AsyncClient() as client:
        auth = AuthToken(client=client)
        res = await client.patch(
            f"{settings.BASE_URL}/project/{project_id}/{dataset_id}",
            json=updates_not_none,
            headers=await auth.header(),
        )

        print_json(data=res.json())


# TASK GROUP


@cli.group()
async def task():
    pass


@task.command(name="list")
async def get_task():

    async with httpx.AsyncClient() as client:
        auth = AuthToken(client=client)
        res = await client.get(
            f"{settings.BASE_URL}/task/",
            headers=await auth.header(),
        )

        print_json(data=res.json())


@task.command(name="new")
@click.argument("name", required=True, nargs=1, type=str)
@click.argument(
    "resource_type",
    required=True,
    nargs=1,
)
@click.argument(
    "input_type",
    required=True,
    nargs=1,
)
@click.argument(
    "output_type",
    required=True,
    nargs=1,
)
@click.option(
    "--module",
    nargs=1,
    help=("default args"),
)
@click.option(
    "--default_args",
    nargs=1,
    help=("default args"),
)
@click.option(
    "--subtask_list",
    nargs=1,
    help=("subtask list of the current task"),
)
async def new_task(
    name: str,
    resource_type: str,
    input_type: str,
    output_type: str,
    default_args: Dict = None,
    module: str = "",
    subtask_list: List = None,
):

    # Check that there is no other dataset with the same name
    from fractal.common.models import TaskRead

    async with httpx.AsyncClient() as client:
        auth = AuthToken(client=client)
        res = await client.get(
            f"{settings.BASE_URL}/task/",
            headers=await auth.header(),
        )
        data = res.json()
        task_list = [TaskRead(**item) for item in data]
        existing_task_names = [t.name for t in task_list]

        if name in existing_task_names:
            raise Exception(f"Task name {name} already in use.")

    from fractal.common.models import TaskCreate

    if not default_args:
        default_args = {}
    if not subtask_list:
        subtask_list = []

    resource_type = resource_type.replace("_", " ")
    task = TaskCreate(
        name=name,
        resource_type=resource_type,
        input_type=input_type,
        output_type=output_type,
        default_args=default_args,
        module=module,
        subtask_list=subtask_list,
    )

    async with httpx.AsyncClient() as client:
        auth = AuthToken(client=client)
        res = await client.post(
            f"{settings.BASE_URL}/task/",
            json=task.dict(),
            headers=await auth.header(),
        )

        print_json(data=res.json())


@task.command(name="modify-task")
@click.argument("task_name", required=True, type=str, nargs=1)
@click.option(
    "--new_task_name",
    type=str,
    nargs=1,
)
@click.option(
    "--default_args",
    nargs=1,
)
async def modify_task(
    task_name: str,
    new_task_name: str = None,
    default_args: str = None,
):

    if default_args is None:
        da = None
    else:
        with open(default_args, "r", encoding="utf-8") as m:
            da = json.load(m)

    from fractal.common.models import TaskRead

    # Extract task list
    async with httpx.AsyncClient() as client:
        auth = AuthToken(client=client)
        res = await client.get(
            f"{settings.BASE_URL}/task/",
            headers=await auth.header(),
        )
        data = res.json()
        task_list = [TaskRead(**item) for item in data]

    if new_task_name is not None:
        # Check that there is no other dataset with the same name
        existing_task_names = [t.name for t in task_list]
        if new_task_name in existing_task_names:
            raise Exception(f"Task name {new_task_name} already in use.")

    updates = dict(
        name=new_task_name,
        default_args=da,
    )
    updates_not_none = {
        key: value for key, value in updates.items() if value is not None
    }

    # Find task id
    try:
        task = [t for t in task_list if t.name == task_name][0]
        task_id = task.id
    except IndexError as e:
        raise IndexError(f"Task {task_name} not found", str(e))

    async with httpx.AsyncClient() as client:
        auth = AuthToken(client=client)
        res = await client.patch(
            f"{settings.BASE_URL}/task/{task_id}/",
            json=updates_not_none,
            headers=await auth.header(),
        )

        print_json(data=res.json())


@task.command(name="add-subtask")
@click.argument("parent_task_name", type=str, required=True, nargs=1)
@click.argument("subtask_name", type=str, required=True, nargs=1)
@click.option("--args_json", type=str, nargs=1)
async def add_subtask(
    parent_task_name: str, subtask_name: str, args_json: str = None
):

    # Extract subtask_id
    from fractal.common.models import TaskRead

    async with httpx.AsyncClient() as client:
        auth = AuthToken(client=client)
        res = await client.get(
            f"{settings.BASE_URL}/task/",
            headers=await auth.header(),
        )
        data = res.json()
        task_list = [TaskRead(**item) for item in data]
        parent_task_id = [
            t.id for t in task_list if t.name == parent_task_name
        ][0]
        subtask_id = [t.id for t in task_list if t.name == subtask_name][0]

    # Read args
    if args_json is None:
        args = {}
    else:
        with open(args_json, "r", encoding="utf-8") as json_file:
            args = json.load(json_file)

    # Create subtask
    logging.debug(subtask_id, args)
    subtask = SubtaskCreate(subtask_id=subtask_id, args=args)

    async with httpx.AsyncClient() as client:
        auth = AuthToken(client=client)
        res = await client.post(
            f"{settings.BASE_URL}/task/{parent_task_id}/subtask/",
            json=subtask.dict(),
            headers=await auth.header(),
        )

        print_json(data=res.json())


# APPLY GROUP


@cli.group()
def workflow():
    pass


@workflow.command(name="apply")
@click.argument("project_name", required=True, nargs=1, type=str)
@click.argument(
    "input_dataset_name",
    required=True,
    type=str,
    nargs=1,
)
@click.option(
    "--output_dataset_name",
    nargs=1,
    type=str,
    help=("output dataset name"),
)
@click.argument("workflow_name", required=True, nargs=1, type=str)
@click.option(
    "--overwrite_input",
    default=False,
    nargs=1,
)
async def apply_workflow(
    project_name: str,
    input_dataset_name: str,
    output_dataset_name: str,
    workflow_name: str,
    overwrite_input: bool,
):

    from fractal.common.models import TaskRead

    # Extract IDs for project and two datasets
    async with httpx.AsyncClient() as client:
        auth = AuthToken(client=client)
        res = await client.get(
            f"{settings.BASE_URL}/project/",
            headers=await auth.header(),
        )
        projects = res.json()
        # Find project
        try:
            project = [p for p in projects if p["name"] == project_name][0]
            project_id = project["id"]
        except IndexError as e:
            raise IndexError(f"Project {project_name} not found", str(e))
        # Get dataset list
        dataset_list = [
            project["dataset_list"]
            for project in projects
            if project["id"] == project_id
        ][0]
        # Find I/O dataset
        try:
            input_dataset = [
                ds for ds in dataset_list if ds["name"] == input_dataset_name
            ][0]
            input_dataset_id = input_dataset["id"]
        except IndexError as e:
            raise IndexError(f"Dataset {input_dataset_name} not found", str(e))
        try:
            output_dataset = [
                ds for ds in dataset_list if ds["name"] == output_dataset_name
            ][0]
            output_dataset_id = output_dataset["id"]
        except IndexError as e:
            raise IndexError(
                f"Dataset {output_dataset_name} not found", str(e)
            )
    # Extract ID for task
    async with httpx.AsyncClient() as client:
        auth = AuthToken(client=client)
        res = await client.get(
            f"{settings.BASE_URL}/task/",
            headers=await auth.header(),
        )
        data = res.json()
        task_list = [TaskRead(**item) for item in data]
        # Find task id
        try:
            workflow = [t for t in task_list if t.name == workflow_name][0]
            workflow_id = workflow.id
        except IndexError as e:
            raise IndexError(f"Task {workflow_name} not found", str(e))

    # Apply workflow
    async with httpx.AsyncClient() as client:
        auth = AuthToken(client=client)

        from fractal.common.models import ApplyWorkflow

        workflow = ApplyWorkflow(
            project_id=project_id,
            input_dataset_id=input_dataset_id,
            output_dataset_id=output_dataset_id,
            workflow_id=workflow_id,
            overwrite_input=overwrite_input,
        )

        res = await client.post(
            f"{settings.BASE_URL}/project/apply/",
            json=workflow.dict(),
            headers=await auth.header(),
        )

        logging.debug(res.json())
