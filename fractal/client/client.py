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
from typing import Any
from typing import Dict
from typing import Optional
import asyncclick as click
import httpx
from rich import print_json
from rich.table import Table
from rich.console import Console
from devtools import debug  # FIXME remove noqa

from ._auth import AuthToken
from .config import settings

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
        #debug(res.json())
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
        table.add_column("ID", style="cyan", no_wrap=True)
        table.add_column("Name", style="magenta")
        table.add_column("Proj. Dir.", justify="right", style="green")
        table.add_column("Dataset list", style="white")
        table.add_column("Read only", justify="center")

        for p in project_list:
            if p.read_only==True:
                p_read_only = "✅"
            p_read_only = "❌"

            p_dataset_list = str([dataset.name for dataset in p.dataset_list])

            table.add_row(str(p.id), p.name, p.project_dir, 
                          str(p_dataset_list),  p_read_only)
        
        console.print(table)


@project.command(name="add-dataset")
@click.argument("project_id", required=True, nargs=1)
@click.argument(
    "name",
    required=True,
    nargs=1,
)
@click.argument(
    "meta",
    required=True,
    nargs=1,
)
@click.option(
    "--type",
    required=True,
    nargs=1,
    default="default",
    help=(
        "The type of objects into the dataset"
    ),
)
@click.option(
    "--read_only",
    required=True,
    nargs=1,
    default="default",
    help=(
        "Behaviour of the dataset"
    ),
)

async def add_dataset(project_id: str,
                      name:str,
                      type: Optional[str],
                      meta: Dict[str,Any] = None,
                      read_only: Optional[bool] = False,
 ) -> None:
    """
    Add an existing dataset to an exisisting project

    project_id (int): project id

    dataset_id (str): dataset id
    """
    from fractal.common.models import DatasetBase

    dataset = DatasetBase(name=name, project_id=project_id,
                         type=type, meta=meta, read_only=read_only)
    async with httpx.AsyncClient() as client:
        auth = AuthToken(client=client)
        res = await client.post(
            f"{settings.BASE_URL}/{project_id}/",
            json=dataset.dict(),
            headers=await auth.header(),
        )
        debug(res.json())

