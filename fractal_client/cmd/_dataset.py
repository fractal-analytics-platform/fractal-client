import json
from typing import Optional

from rich.table import Table

from ..authclient import AuthClient
from ..config import settings
from ..interface import BaseInterface
from ..interface import PrintInterface
from ..interface import RichConsoleInterface
from ..interface import RichJsonInterface
from ..response import check_response


def post_dataset(
    client: AuthClient,
    *,
    project_id: int,
    dataset_name: str,
    metadata: Optional[str] = None,
    type: Optional[str] = None,
    batch: bool = False,
    make_read_only: bool = False,
) -> RichJsonInterface:
    """
    Arguments:
        project_id: ID of project to add the new dataset to
        dataset_name: Name of new dataset
        metadata: Path to file containing dataset metadata in JSON format.
        type: Dataset type.
        make_read_only: Make the new dataset read-only.
    """
    if metadata is None:
        meta = {}
    else:
        with open(metadata, "r") as f:
            meta = json.load(f)

    dataset = dict(name=dataset_name, meta=meta, read_only=make_read_only)
    if type:
        dataset["type"] = type

    res = client.post(
        f"{settings.BASE_URL}/project/{project_id}/dataset/",
        json=dataset,
    )
    new_dataset = check_response(res, expected_status_code=201)
    if batch:
        return PrintInterface(retcode=0, data=new_dataset["id"])
    else:
        return RichJsonInterface(retcode=0, data=new_dataset)


def post_resource(
    client: AuthClient,
    *,
    project_id: int,
    dataset_id: int,
    path: str,
    batch: bool = False,
) -> BaseInterface:

    res = client.post(
        (
            f"{settings.BASE_URL}/project/{project_id}/"
            f"dataset/{dataset_id}/resource/"
        ),
        json=dict(path=path),
    )
    new_resource = check_response(res, expected_status_code=201)
    if batch:
        return PrintInterface(retcode=0, data=new_resource["id"])
    else:
        return RichJsonInterface(retcode=0, data=new_resource)


def delete_resource(
    client: AuthClient,
    *,
    project_id: int,
    dataset_id: int,
    resource_id: int,
) -> BaseInterface:
    res = client.delete(
        (
            f"{settings.BASE_URL}/project/{project_id}/"
            f"dataset/{dataset_id}/resource/{resource_id}/"
        )
    )
    check_response(res, expected_status_code=204)
    return PrintInterface(retcode=0, data="")


def patch_dataset(
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
    dataset_update = {}
    if new_name:
        dataset_update["name"] = new_name
    if new_type:
        dataset_update["type"] = new_type
    if meta_file:
        with open(meta_file, "r") as f:
            meta = json.load(f)
        dataset_update.update(meta=meta)
    if make_read_only:
        dataset_update["read_only"] = True
    if remove_read_only:
        dataset_update["read_only"] = False

    if not dataset_update:
        return PrintInterface(retcode=1, data="Nothing to update")

    res = client.patch(
        (
            f"{settings.BASE_URL}/project/{project_id}/"
            f"dataset/{dataset_id}/"
        ),
        json=dataset_update,
    )
    data = check_response(res, expected_status_code=200)
    return RichJsonInterface(retcode=0, data=data)


def get_dataset(
    client: AuthClient, *, project_id: int, dataset_id: int
) -> BaseInterface:
    res = client.get(
        f"{settings.BASE_URL}/project/{project_id}/dataset/{dataset_id}/"
    )
    from rich.console import Group

    dataset = check_response(res, expected_status_code=200)

    table = Table(title="Dataset")
    table.add_column("ID", style="cyan", no_wrap=True)
    table.add_column("Name", justify="right", style="green")
    table.add_column("Type", style="white")
    table.add_column("Meta", justify="center")
    table.add_column("Read only", justify="center")

    table.add_row(
        str(dataset["id"]),
        dataset["name"],
        dataset["type"],
        json.dumps(dataset["meta"], indent=2),
        "✅" if dataset["read_only"] else "❌",
    )
    table_res = Table(title="Resources")
    table_res.add_column("Path", justify="center", style="yellow")
    table_res.add_column("ID", justify="center", style="yellow")
    table_res.add_column("Dataset ID", justify="center", style="yellow")
    for r in dataset["resource_list"]:
        table_res.add_row(r["path"], str(r["id"]), str(r["dataset_id"]))
    group = Group(table, table_res)
    return RichConsoleInterface(retcode=0, data=group)


def delete_dataset(
    client: AuthClient, *, project_id: int, dataset_id: int
) -> PrintInterface:

    res = client.delete(
        f"{settings.BASE_URL}/project/{project_id}/dataset/{dataset_id}/"
    )
    check_response(res, expected_status_code=204)
    return PrintInterface(retcode=0, data="")


def get_dataset_history(
    client: AuthClient, *, project_id: int, dataset_id: int
) -> BaseInterface:
    res = client.get(
        f"{settings.BASE_URL}/project/{project_id}/dataset/{dataset_id}/"
        "export_history/"
    )
    history_workflow = check_response(res, expected_status_code=200)
    return RichJsonInterface(retcode=0, data=history_workflow)


def get_dataset_status(
    client: AuthClient, *, project_id: int, dataset_id: int
) -> BaseInterface:
    res = client.get(
        f"{settings.BASE_URL}/project/{project_id}/dataset/{dataset_id}/"
        "status/"
    )
    dataset_status = check_response(res, expected_status_code=200)
    return RichJsonInterface(retcode=0, data=dataset_status)
