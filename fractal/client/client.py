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
import asyncclick as click
import httpx
from devtools import debug  # FIXME remove noqa

from ._auth import AuthToken
from .config import settings


@click.group()
async def cli():
    pass


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
    from ..server.app.models import Project

    project = Project(name=name, project_dir=path)
    async with httpx.AsyncClient() as client:
        auth = AuthToken(client=client)
        res = await client.post(
            f"{settings.BASE_URL}/project/",
            json=project.dict(),
            headers=await auth.header(),
        )
        debug(res.json())


@project.command(name="list")
async def project_list():
    async with httpx.AsyncClient() as client:
        auth = AuthToken(client=client)
        res = await client.get(
            f"{settings.BASE_URL}/project/",
            headers=await auth.header(),
        )
        debug(res.json())


@cli.command(name="login")
async def login():
    async with httpx.AsyncClient() as client:
        auth = AuthToken(client=client)
        await auth()
        debug(await auth.header())
