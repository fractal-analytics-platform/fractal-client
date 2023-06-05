import json
import os
from typing import Optional

from rich.table import Table

from ..authclient import AuthClient
from ..common.schemas import DatasetCreate
from ..common.schemas import DatasetRead
from ..common.schemas import DatasetUpdate
from ..common.schemas import ResourceCreate
from ..common.schemas import ResourceRead
from ..config import settings
from ..interface import BaseInterface
from ..interface import PrintInterface
from ..interface import RichConsoleInterface
from ..interface import RichJsonInterface
from ..response import check_response


async def post_dataset(
    client: AuthClient,
    *,
    project_id: int,
    dataset_name: str,
    metadata: Optional[str] = None,
    type: Optional[str] = None,
    batch: bool = False,
) -> RichJsonInterface:
    """
    Arguments:
        project_id: ID of project to add the new dataset to
        dataset_name: Name of new dataset
        metadata: Path to file containing dataset metadata in JSON format.
        type: Dataset type.
    """
    if metadata is None:
        meta = {}
    else:
        with open(metadata, "r") as f:
            meta = json.load(f)

    dataset_dict = dict(name=dataset_name, meta=meta)
    if type:
        dataset_dict["type"] = type
    dataset = DatasetCreate(**dataset_dict)

    res = await client.post(
        f"{settings.BASE_URL}/project/{project_id}/dataset/",
        json=dataset.dict(exclude_unset=True),
    )
    new_dataset = check_response(
        res, expected_status_code=201, coerce=DatasetRead
    )
    if batch:
        return PrintInterface(retcode=0, data=new_dataset.id)
    else:
        return RichJsonInterface(retcode=0, data=new_dataset.dict())


async def post_resource(
    client: AuthClient,
    *,
    project_id: int,
    dataset_id: int,
    path: str,
    batch: bool = False,
) -> BaseInterface:

    # Check that path is absolute, which is needed for when the server submits
    # tasks as a different user
    if not os.path.isabs(path):
        msg = f"{path=} is not an absolute path"
        raise ValueError(msg)

    resource = ResourceCreate(path=path)

    res = await client.post(
        (
            f"{settings.BASE_URL}/project/{project_id}/"
            f"dataset/{dataset_id}/resource/"
        ),
        json=resource.dict(),
    )
    new_resource = check_response(
        res, expected_status_code=201, coerce=ResourceRead
    )
    if batch:
        return PrintInterface(retcode=0, data=new_resource.id)
    else:
        return RichJsonInterface(retcode=0, data=new_resource.dict())


async def delete_resource(
    client: AuthClient,
    *,
    project_id: int,
    dataset_id: int,
    resource_id: int,
) -> BaseInterface:
    res = await client.delete(
        (
            f"{settings.BASE_URL}/project/{project_id}/"
            f"dataset/{dataset_id}/resource/{resource_id}"
        )
    )
    check_response(res, expected_status_code=204)
    return PrintInterface(retcode=0, data="")


async def patch_dataset(
    client: AuthClient,
    *,
    project_id: int,
    dataset_id: int,
    new_name: Optional[str] = None,
    new_type: Optional[str] = None,
    meta_file: Optional[str] = None,
    make_read_only: bool = False,
    remove_read_only: bool = False,
) -> BaseInterface:

    # Prepare payload
    dataset_update_dict = {}
    if new_name:
        dataset_update_dict["name"] = new_name
    if new_type:
        dataset_update_dict["type"] = new_type
    if meta_file:
        with open(meta_file, "r") as f:
            meta = json.load(f)
        dataset_update_dict.update(meta=meta)
    if make_read_only:
        dataset_update_dict["read_only"] = True
    if remove_read_only:
        dataset_update_dict["read_only"] = False
    dataset_update = DatasetUpdate(**dataset_update_dict)
    payload = dataset_update.dict(exclude_unset=True)

    if not payload:
        return PrintInterface(retcode=1, data="Nothing to update")

    res = await client.patch(
        (f"{settings.BASE_URL}/project/{project_id}/" f"dataset/{dataset_id}"),
        json=payload,
    )
    new_dataset = check_response(
        res, expected_status_code=200, coerce=DatasetRead
    )
    return RichJsonInterface(retcode=0, data=new_dataset.dict())


async def get_dataset(
    client: AuthClient, *, project_id: int, dataset_id: int
) -> BaseInterface:
    res = await client.get(
        f"{settings.BASE_URL}/project/{project_id}/dataset/{dataset_id}"
    )
    from rich.console import Group

    dataset = check_response(res, expected_status_code=200, coerce=DatasetRead)

    table = Table(title="Dataset")
    table.add_column("ID", style="cyan", no_wrap=True)
    table.add_column("Name", justify="right", style="green")
    table.add_column("Type", style="white")
    table.add_column("Meta", justify="center")
    table.add_column("Read only", justify="center")

    table.add_row(
        str(dataset.id),
        dataset.name,
        dataset.type,
        json.dumps(dataset.meta, indent=2),
        "✅" if dataset.read_only else "❌",
    )
    table_res = Table(title="Resources")
    table_res.add_column("Path", justify="center", style="yellow")
    table_res.add_column("ID", justify="center", style="yellow")
    table_res.add_column("Dataset ID", justify="center", style="yellow")
    for r in dataset.resource_list:
        table_res.add_row(r.path, str(r.id), str(r.dataset_id))
    group = Group(table, table_res)
    return RichConsoleInterface(retcode=0, data=group)


async def delete_dataset(
    client: AuthClient, *, project_id: int, dataset_id: int
) -> PrintInterface:

    res = await client.delete(
        f"{settings.BASE_URL}/project/{project_id}/dataset/{dataset_id}"
    )
    check_response(res, expected_status_code=204)
    return PrintInterface(retcode=0, data="")
