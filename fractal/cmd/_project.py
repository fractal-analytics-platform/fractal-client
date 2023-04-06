import json
import logging
from typing import Optional
from typing import Union

from rich.table import Table

from ..authclient import AuthClient
from ..common.schemas import DatasetCreate
from ..common.schemas import DatasetRead
from ..common.schemas import ProjectCreate
from ..common.schemas import ProjectRead
from ..common.schemas import ProjectUpdate
from ..config import settings
from ..interface import BaseInterface
from ..interface import PrintInterface
from ..interface import RichConsoleInterface
from ..interface import RichJsonInterface
from ..response import check_response


async def project_create(
    client: AuthClient,
    name: str,
    dataset: Optional[str] = None,
    batch: bool = False,
    **kwargs,
) -> BaseInterface:
    # Prepare a ProjectCreate request body
    project_dict = dict(name=name)
    if dataset:
        project_dict["default_dataset_name"] = dataset
    project = ProjectCreate(**project_dict)
    logging.info(project)
    # Send API request
    res = await client.post(
        f"{settings.BASE_URL}/project/",
        json=project.dict(exclude_unset=True),
    )
    project = check_response(res, expected_status_code=201, coerce=ProjectRead)
    if batch:
        if len(project.dataset_list) > 1:
            msg = (
                f"Created project with {len(project.dataset_list)}>1 "
                "datasets, cannot use --batch to provide standard output."
            )
            raise ValueError(msg)
        dataset_id = project.dataset_list[0].id
        return PrintInterface(retcode=0, data=f"{project.id} {dataset_id}")
    else:
        return RichJsonInterface(retcode=0, data=project.dict())


async def project_list(client: AuthClient, **kwargs) -> RichConsoleInterface:

    res = await client.get(f"{settings.BASE_URL}/project/")
    res = check_response(res, expected_status_code=200)

    projects = [ProjectRead(**item) for item in res]

    table = Table(title="Project List")
    table.add_column("ID", style="cyan", no_wrap=True)
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
            str(p_dataset_list),
            read_only_icon,
        )

    return RichConsoleInterface(retcode=0, data=table)


async def project_show(
    client: AuthClient, project_id: int, **kwargs
) -> RichJsonInterface:
    res = await client.get(
        f"{settings.BASE_URL}/project/{project_id}",
    )
    project = check_response(res, expected_status_code=200, coerce=ProjectRead)
    return RichJsonInterface(retcode=0, data=project.dict())


async def project_add_dataset(
    client: AuthClient,
    project_id: int,
    dataset_name: str,
    metadata_filename: Optional[str] = None,
    type: Optional[str] = None,
    batch: bool = False,
    **kwargs,
) -> RichJsonInterface:

    if metadata_filename is None:
        meta = {}
    else:
        meta = json.loads(metadata_filename)

    dataset_dict = dict(name=dataset_name, meta=meta)
    if type:
        dataset_dict["type"] = type
    dataset = DatasetCreate(**dataset_dict)

    res = await client.post(
        f"{settings.BASE_URL}/project/{project_id}/",
        json=dataset.dict(exclude_unset=True),
    )
    new_dataset = check_response(
        res, expected_status_code=201, coerce=DatasetRead
    )
    if batch:
        return PrintInterface(retcode=0, data=new_dataset.id)
    else:
        return RichJsonInterface(retcode=0, data=new_dataset.dict())


async def project_delete(
    client: AuthClient, project_id: int, **kwargs
) -> PrintInterface:

    res = await client.delete(
        f"{settings.BASE_URL}/project/{project_id}",
    )
    check_response(res, expected_status_code=204)
    return PrintInterface(retcode=0, data="")


async def project_edit(
    client: AuthClient,
    project_id: int,
    new_name: Optional[str] = None,
    make_read_only: bool = False,
    remove_read_only: bool = False,
    **kwargs,
) -> Union[RichJsonInterface, PrintInterface]:
    update = {}
    if new_name:
        update["name"] = new_name
    if make_read_only:
        update["read_only"] = True
    if remove_read_only:
        update["read_only"] = False

    project_update = ProjectUpdate(**update)  # validation
    payload = project_update.dict(exclude_unset=True)

    if not payload:
        return PrintInterface(retcode=1, data="Nothing to update")

    res = await client.patch(
        f"{settings.BASE_URL}/project/{project_id}", json=payload
    )
    new_project = check_response(
        res, expected_status_code=200, coerce=ProjectRead
    )
    return RichJsonInterface(retcode=0, data=new_project.dict())
