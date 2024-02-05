import json
from typing import Optional

from ..authclient import AuthClient
from ..config import settings
from ..interface import Interface
from ..response import check_response
from ._aux_trim_output import _simplify_dataset
from ._aux_trim_output import _simplify_resource


def post_dataset(
    client: AuthClient,
    *,
    project_id: int,
    dataset_name: str,
    metadata: Optional[str] = None,
    type: Optional[str] = None,
    batch: bool = False,
    verbose: bool = False,
    make_read_only: bool = False,
) -> Interface:
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
        return Interface(retcode=0, data=new_dataset["id"])
    elif verbose:
        return Interface(retcode=0, data=new_dataset)
    else:
        return Interface(retcode=0, data=_simplify_dataset(new_dataset))


def post_resource(
    client: AuthClient,
    *,
    project_id: int,
    dataset_id: int,
    path: str,
    batch: bool = False,
    verbose: bool = False,
) -> Interface:

    res = client.post(
        (
            f"{settings.BASE_URL}/project/{project_id}/"
            f"dataset/{dataset_id}/resource/"
        ),
        json=dict(path=path),
    )
    new_resource = check_response(res, expected_status_code=201)
    if batch:
        return Interface(retcode=0, data=new_resource["id"])
    elif verbose:
        return Interface(retcode=0, data=new_resource)
    else:
        return Interface(retcode=0, data=_simplify_resource(new_resource))


def delete_resource(
    client: AuthClient,
    *,
    project_id: int,
    dataset_id: int,
    resource_id: int,
) -> Interface:
    res = client.delete(
        (
            f"{settings.BASE_URL}/project/{project_id}/"
            f"dataset/{dataset_id}/resource/{resource_id}/"
        )
    )
    check_response(res, expected_status_code=204)
    return Interface(retcode=0, data="")


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
    batch: bool = False,
    verbose: bool = False,
) -> Interface:

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
        return Interface(retcode=1, data="Nothing to update")

    res = client.patch(
        (
            f"{settings.BASE_URL}/project/{project_id}/"
            f"dataset/{dataset_id}/"
        ),
        json=dataset_update,
    )
    data = check_response(res, expected_status_code=200)
    if batch:
        return Interface(retcode=0, data=data["id"])
    elif verbose:
        return Interface(retcode=0, data=data)
    else:
        return Interface(retcode=0, data=_simplify_dataset(data))


def get_dataset(
    client: AuthClient,
    *,
    project_id: int,
    dataset_id: int,
    batch: bool = False,
    verbose: bool = False,
) -> Interface:
    res = client.get(
        f"{settings.BASE_URL}/project/{project_id}/dataset/{dataset_id}/"
    )
    dataset = check_response(res, expected_status_code=200)
    if batch:
        return Interface(retcode=0, data=dataset["id"])
    elif verbose:
        return Interface(retcode=0, data=dataset)
    else:
        return Interface(retcode=0, data=_simplify_dataset(dataset))


def delete_dataset(
    client: AuthClient, *, project_id: int, dataset_id: int
) -> Interface:

    res = client.delete(
        f"{settings.BASE_URL}/project/{project_id}/dataset/{dataset_id}/"
    )
    check_response(res, expected_status_code=204)
    return Interface(retcode=0, data="")


def get_dataset_history(
    client: AuthClient,
    *,
    project_id: int,
    dataset_id: int,
    batch: bool = False,  # FIXME
    verbose: bool = False,  # FIXME
) -> Interface:
    res = client.get(
        f"{settings.BASE_URL}/project/{project_id}/dataset/{dataset_id}/"
        "export_history/"
    )
    history_workflow = check_response(res, expected_status_code=200)
    return Interface(retcode=0, data=history_workflow)


def get_dataset_status(
    client: AuthClient, *, project_id: int, dataset_id: int
) -> Interface:
    res = client.get(
        f"{settings.BASE_URL}/project/{project_id}/dataset/{dataset_id}/"
        "status/"
    )
    dataset_status = check_response(res, expected_status_code=200)
    return Interface(retcode=0, data=dataset_status)
