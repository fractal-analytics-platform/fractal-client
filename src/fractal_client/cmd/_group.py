from ..authclient import AuthClient
from ..interface import Interface
from ..response import check_response


def group_list(
    client: AuthClient,
    *,
    user_ids: bool = False,
    batch: bool = False,
):
    query_params = "?user_ids=true" if user_ids else ""
    res = client.get(f"auth/group/{query_params}")
    data = check_response(res, expected_status_code=200)
    if batch:
        return Interface(
            retcode=0, data=" ".join([str(d["id"]) for d in data])
        )
    else:
        return Interface(retcode=0, data=data)


def group_get(client: AuthClient, *, group_id: int):
    res = client.get(f"auth/group/{group_id}/")
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
        "auth/group/",
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
    new_viewer_paths: list[str],
):
    res = client.patch(
        f"auth/group/{group_id}/",
        json=dict(viewer_paths=new_viewer_paths),
    )
    data = check_response(res, expected_status_code=200)
    return Interface(retcode=0, data=data)


def group_add_user(
    client: AuthClient,
    *,
    group_id: int,
    user_id: int,
):
    res = client.post(f"auth/group/{group_id}/add-user/{user_id}/")
    data = check_response(res, expected_status_code=200)
    return Interface(retcode=0, data=data)


def group_remove_user(
    client: AuthClient,
    *,
    group_id: int,
    user_id: int,
):
    res = client.post(f"auth/group/{group_id}/remove-user/{user_id}/")
    data = check_response(res, expected_status_code=200)
    return Interface(retcode=0, data=data)
