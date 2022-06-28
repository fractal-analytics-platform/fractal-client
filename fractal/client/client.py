import asyncclick as click
import httpx

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
        res = await client.post(
            f"{settings.BASE_URL}/project/", json=project.dict()
        )
        from devtools import debug

        debug(res.json())
