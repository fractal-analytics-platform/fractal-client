from ..authclient import AuthClient
from ..interface import Interface
from ..response import check_response


def post_dataset(
    client: AuthClient,
    *,
    project_id: int,
    dataset_name: str,
    project_dir: str | None = None,
    zarr_subfolder: str | None = None,
    batch: bool = False,
) -> Interface:
    """
    Arguments:
        project_id: ID of project to add the new dataset to
        dataset_name: Name of new dataset
        project_dir:
        zarr_subfolder:
        batch:
    """
    dataset = dict(name=dataset_name)
    if project_dir is not None:
        dataset["project_dir"] = project_dir
    if zarr_subfolder is not None:
        dataset["zarr_subfolder"] = zarr_subfolder

    res = client.post(
        f"api/v2/project/{project_id}/dataset/",
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
    new_name: str,
) -> Interface:
    res = client.patch(
        (f"api/v2/project/{project_id}/" f"dataset/{dataset_id}/"),
        json=dict(name=new_name),
    )
    data = check_response(res, expected_status_code=200)
    return Interface(retcode=0, data=data)


def get_dataset(
    client: AuthClient, *, project_id: int, dataset_id: int
) -> Interface:
    res = client.get(f"api/v2/project/{project_id}/dataset/{dataset_id}/")
    return Interface(retcode=0, data=res.json())


def delete_dataset(
    client: AuthClient, *, project_id: int, dataset_id: int
) -> Interface:
    res = client.delete(f"api/v2/project/{project_id}/dataset/{dataset_id}/")
    check_response(res, expected_status_code=204)
    return Interface(retcode=0, data="")
