from ..authclient import AuthClient
from ..interface import Interface
from ..response import check_response


def post_project(
    client: AuthClient,
    *,
    name: str,
    batch: bool = False,
) -> Interface:
    # Prepare a ProjectCreate request body
    payload = dict(name=name)
    # Send API request
    res = client.post("api/v2/project/", json=payload)
    project = check_response(res, expected_status_code=201)
    if batch:
        return Interface(retcode=0, data=project["id"])
    else:
        return Interface(retcode=0, data=project)


def get_project_list(client: AuthClient) -> Interface:

    res = client.get("api/v2/project/")
    projects = check_response(res, expected_status_code=200)
    return Interface(retcode=0, data=projects)


def get_project(client: AuthClient, *, project_id: int) -> Interface:
    res = client.get(f"api/v2/project/{project_id}/")
    project = check_response(res, expected_status_code=200)
    return Interface(retcode=0, data=project)


def delete_project(client: AuthClient, *, project_id: int) -> Interface:

    res = client.delete(f"api/v2/project/{project_id}/")
    check_response(res, expected_status_code=204)
    return Interface(retcode=0, data="")


def patch_project(
    client: AuthClient,
    *,
    project_id: int,
    new_name: str | None = None,
) -> Interface:
    project_update = {}
    if new_name:
        project_update["name"] = new_name

    res = client.patch(f"api/v2/project/{project_id}/", json=project_update)
    new_project = check_response(res, expected_status_code=200)
    return Interface(retcode=0, data=new_project)
