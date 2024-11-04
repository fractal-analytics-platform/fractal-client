from ..authclient import AuthClient
from ..config import settings
from ..interface import Interface
from ..response import check_response


def group_list(
    client: AuthClient, *, user_ids: bool = False, batch: bool = False
):
    query_params = "?user_ids=true" if user_ids else ""
    res = client.get(f"{settings.FRACTAL_SERVER}/auth/group/{query_params}")
    data = check_response(res, expected_status_code=200)
    if batch:
        return Interface(
            retcode=0, data=" ".join([str(d["id"]) for d in data])
        )
    else:
        return Interface(retcode=0, data=data)


def group_get(client: AuthClient, *, group_id: int):
    res = client.get(f"{settings.FRACTAL_SERVER}/auth/group/{group_id}/")
    data = check_response(res, expected_status_code=200)
    return Interface(retcode=0, data=data)


def group_new(
    client: AuthClient,
    *,
    name: str,
    viewer_paths: list[str] | None = None,
    batch: bool = False,
):
    request_body = dict(name=name)
    if viewer_paths is not None:
        request_body["viewer_paths"] = viewer_paths

    res = client.post(
        f"{settings.FRACTAL_SERVER}/auth/group/",
        json=request_body,
    )
    data = check_response(res, expected_status_code=201)
    if batch:
        return Interface(retcode=0, data=data["id"])
    else:
        return Interface(retcode=0, data=data)


def group_update(
    client: AuthClient,
    *,
    group_id: int,
    new_user_ids: list[int] | None = None,
    new_viewer_paths: list[str] | None = None,
):

    request_body = dict()
    if new_viewer_paths is not None:
        request_body["viewer_paths"] = new_viewer_paths
    if new_user_ids is not None:
        request_body["new_user_ids"] = new_user_ids

    res = client.patch(
        f"{settings.FRACTAL_SERVER}/auth/group/{group_id}/",
        json=request_body,
    )
    data = check_response(res, expected_status_code=200)
    return Interface(retcode=0, data=data)
