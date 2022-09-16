import json
import logging
from typing import Optional

from rich import print_json
from rich.console import Console
from rich.table import Table

from ...common.models import DatasetCreate
from ...common.models import DatasetRead
from ...common.models import ProjectCreate
from ...common.models import ProjectRead
from ..authclient import AuthClient
from ..config import settings
from ..response import check_response


async def project_create(
    client: AuthClient,
    name: str,
    path: str,
    dataset: Optional[str] = None,
    batch: bool = False,
    **kwargs,
) -> None:
    project = ProjectCreate(
        name=name, project_dir=path, default_dataset_name=dataset
    )
    logging.info(project)
    res = await client.post(
        f"{settings.BASE_URL}/project/",
        json=project.dict(),
    )
    project = check_response(res, expected_status_code=201, coerce=ProjectRead)
    if batch:
        print(project.id)
    else:
        print(project)


async def project_list(client: AuthClient, **kwargs) -> None:
    res = await client.get(
        f"{settings.BASE_URL}/project/",
    )
    data = res.json()
    projects = [ProjectRead(**item) for item in data]

    table = Table(title="Project List")
    table.add_column("Id", style="cyan", no_wrap=True)
    table.add_column("Name", style="magenta")
    table.add_column("Proj. Dir.", justify="right", style="green")
    table.add_column("Dataset list", style="white")
    table.add_column("Read only", justify="center")

    for p in projects:
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

    console = Console()
    console.print(table)


async def project_add_dataset(
    client: AuthClient,
    project_id: int,
    dataset_name: str,
    metadata_filename: Optional[str] = None,
    type: Optional[str] = None,
    **kwargs,
):

    if metadata_filename is None:
        meta = {}
    else:
        meta = json.loads(metadata_filename)

    dataset = DatasetCreate(
        name=dataset_name, project_id=project_id, type=type, meta=meta
    )

    res = await client.post(
        f"{settings.BASE_URL}/project/{project_id}/",
        json=dataset.dict(),
    )
    new_dataset = check_response(
        res, expected_status_code=201, coerce=DatasetRead
    )
    print_json(data=new_dataset.dict())
