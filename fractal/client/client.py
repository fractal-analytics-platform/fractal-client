"""
Copyright 2022 (C) Friedrich Miescher Institute for Biomedical Research and
University of Zurich

Original authors:
Jacopo Nespolo <jacopo.nespolo@exact-lab.it>

This file is part of Fractal and was originally developed by eXact lab S.r.l.
<exact-lab.it> under contract with Liberali Lab from the Friedrich Miescher
Institute for Biomedical Research and Pelkmans Lab from the University of
Zurich.
"""
import json
from typing import Any
from typing import Dict
from typing import Optional

import asyncclick as click
import httpx
from devtools import debug  # FIXME remove noqa
from rich import print_json
from rich.console import Console
from rich.table import Table

from ._auth import AuthToken
from .config import settings
from fractal.common.models import ResourceRead

console = Console()


@click.group()
async def cli():
    pass


@cli.command(name="login")
async def login():
    async with httpx.AsyncClient() as client:
        auth = AuthToken(client=client)
        await auth()
        debug(await auth.header())


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
        # debug(res.json())
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
            if p.read_only:
                p_read_only = "✅"
            p_read_only = "❌"

            p_dataset_list = str([dataset.name for dataset in p.dataset_list])

            table.add_row(
                str(p.id),
                p.name,
                p.project_dir,
                str(p_dataset_list),
                p_read_only,
            )

        console.print(table)


@project.command(name="add-dataset")
@click.argument("project_id", required=True, nargs=1)
@click.argument(
    "name_dataset",
    required=True,
    nargs=1,
)
@click.argument(
    "meta",
    required=True,
    type=click.File("rb"),
    nargs=1,
)
@click.option(
    "--type",
    required=True,
    nargs=1,
    default="zarr",
    help=("The type of objects into the dataset"),
)
async def add_dataset(
    project_id: str,
    name_dataset: str,
    type: Optional[str],
    meta: Dict[str, Any],
) -> None:
    """
    Add an existing dataset to an exisisting project
    """
    from fractal.common.models import DatasetCreate

    meta_json = json.load(meta)

    dataset = DatasetCreate(
        name=name_dataset, project_id=project_id, type=type, meta=meta_json
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
@click.argument("project_id", required=True, type=int, nargs=1)
@click.argument("dataset_name", required=True, type=str, nargs=1)
async def dataset_show(project_id: int, dataset_name: str) -> None:
    """
    Show details of an exisisting dataset
    """

    async with httpx.AsyncClient() as client:
        auth = AuthToken(client=client)
        res = await client.get(
            f"{settings.BASE_URL}/project/",
            headers=await auth.header(),
        )
        projects = res.json()
        dataset_list = [
            project["dataset_list"]
            for project in projects
            if project["id"] == project_id
        ][0]
        dataset = [ds for ds in dataset_list if ds["name"] == dataset_name][0]

        table = Table(title="Dataset")
        table.add_column("Id", style="cyan", no_wrap=True)
        table.add_column("Name", justify="right", style="green")
        table.add_column("Type", style="white")
        table.add_column("Meta", justify="center")

        table.add_row(
            str(dataset["id"]),
            dataset["name"],
            dataset["type"],
            str(dataset["meta"]),
        )

        console.print(table)


@dataset.command(name="add-resource")
@click.argument("project_id", required=True, type=int, nargs=1)
@click.argument("dataset_id", required=True, type=int, nargs=1)
@click.argument("path", required=True, type=str, nargs=1)
@click.option("--glob_pattern", required=True, type=str, default="", nargs=1)
async def add_resource(
    project_id: int, dataset_id: int, path: str, glob_pattern: str
):

    """
    Add a new resource to an exisisting dataset
    """

    from fractal.common.models import ResourceCreate

    resource = ResourceCreate(path=path, glob_pattern=glob_pattern)

    async with httpx.AsyncClient() as client:
        auth = AuthToken(client=client)
        res = await client.post(
            f"{settings.BASE_URL}/project/{project_id}/{dataset_id}",
            json=resource.dict(),
            headers=await auth.header(),
        )

        print_json(data=res.json())


@dataset.command(name="show-resources")
@click.argument("project_id", required=True, type=int, nargs=1)
@click.argument("dataset_id", required=True, type=int, nargs=1)
async def get_resource(
    project_id: int,
    dataset_id: int,
):

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

        for r in resource_list:

            table.add_row(str(r.id), r.path, str(r.dataset_id))

        console.print(table)


@dataset.command(name="modify-dataset")
@click.argument("project_id", required=True, type=int, nargs=1)
@click.argument("dataset_id", required=True, type=int, nargs=1)
@click.option(
    "--name_dataset",
    required=True,
    nargs=1,
)
@click.option(
    "--meta",
    required=True,
    type=click.File("rb"),
    nargs=1,
)
@click.option(
    "--type",
    required=True,
    nargs=1,
    help=("The type of objects into the dataset"),
)
async def modity_update(
    project_id: int,
    dataset_id: int,
    name_dataset: str = "",
    meta: Dict = None,
    type: str = "",
):

    pass
