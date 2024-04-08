import json
from typing import Optional

from ..authclient import AuthClient
from ..config import settings
from ..interface import Interface
from ..response import check_response


def post_dataset(
    client: AuthClient,
    *,
    project_id: int,
    dataset_name: str,
    zarr_dir: str,
    filters: Optional[str] = None,
    batch: bool = False,
) -> Interface:
    """
    Arguments:
        project_id: ID of project to add the new dataset to
        dataset_name: Name of new dataset
        filters: Path to file containing dataset filters in JSON format.
        batch: Dataset filters.
    """
    if filters is None:
        dataset = dict(name=dataset_name, zarr_dir=zarr_dir)
    else:
        with open(filters, "r") as f:
            filters_dict = json.load(f)
        dataset = dict(
            name=dataset_name, filters=filters_dict, zarr_dir=zarr_dir
        )

    res = client.post(
        f"{settings.BASE_URL}/project/{project_id}/dataset/",
        json=dataset,
    )
    new_dataset = check_response(res, expected_status_code=201)
    if batch:
        return Interface(retcode=0, data=new_dataset["id"])
    else:
        return Interface(retcode=0, data=new_dataset)


def patch_dataset(
    client: AuthClient,
    *,
    project_id: int,
    dataset_id: int,
    new_name: Optional[str] = None,
    filters: Optional[str] = None,
) -> Interface:
    # Prepare payload
    dataset_update = {}
    if new_name:
        dataset_update["name"] = new_name
    if filters:
        with open(filters, "r") as f:
            filters_from_file = json.load(f)
        dataset_update["filters"] = filters_from_file

    res = client.patch(
        (
            f"{settings.BASE_URL}/project/{project_id}/"
            f"dataset/{dataset_id}/"
        ),
        json=dataset_update,
    )
    data = check_response(res, expected_status_code=200)
    return Interface(retcode=0, data=data)


def get_dataset(
    client: AuthClient, *, project_id: int, dataset_id: int
) -> Interface:
    res = client.get(
        f"{settings.BASE_URL}/project/{project_id}/dataset/{dataset_id}/"
    )
    return Interface(retcode=0, data=res.json())


def delete_dataset(
    client: AuthClient, *, project_id: int, dataset_id: int
) -> Interface:

    res = client.delete(
        f"{settings.BASE_URL}/project/{project_id}/dataset/{dataset_id}/"
    )
    check_response(res, expected_status_code=204)
    return Interface(retcode=0, data="")
