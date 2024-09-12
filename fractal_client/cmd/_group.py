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


def group_new(client: AuthClient, *, name: str, batch: bool = False):
    res = client.post(
        f"{settings.FRACTAL_SERVER}/auth/group/", json=dict(name=name)
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
    new_user_ids: list[int],
):
    res = client.patch(
        f"{settings.FRACTAL_SERVER}/auth/group/{group_id}/",
        json=dict(new_user_ids=new_user_ids),
    )
    data = check_response(res, expected_status_code=200)
    return Interface(retcode=0, data=data)
