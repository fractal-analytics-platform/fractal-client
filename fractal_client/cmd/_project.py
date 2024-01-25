from typing import Optional
from typing import Union

from rich.table import Table

from ..authclient import AuthClient
from ..config import settings
from ..interface import BaseInterface
from ..interface import PrintInterface
from ..interface import RichConsoleInterface
from ..interface import RichJsonInterface
from ..response import check_response


def post_project(
    client: AuthClient,
    *,
    name: str,
    batch: bool = False,
) -> BaseInterface:
    # Prepare a ProjectCreate request body
    payload = dict(name=name)
    # Send API request
    res = client.post(f"{settings.BASE_URL}/project/", json=payload)
    project = check_response(res, expected_status_code=201)
    if batch:
        return PrintInterface(retcode=0, data=project["id"])
    else:
        return RichJsonInterface(retcode=0, data=project)


def get_project_list(client: AuthClient) -> RichConsoleInterface:

    res = client.get(f"{settings.BASE_URL}/project/")
    projects = check_response(res, expected_status_code=200)

    table = Table(title="Project List")
    table.add_column("ID", style="cyan", no_wrap=True)
    table.add_column("Name", style="magenta")
    table.add_column("Read only", justify="center")

    for p in projects:
        # Map p.read_only (True/False) to read_only_icon (✅/❌)
        if p["read_only"]:
            read_only_icon = "✅"
        else:
            read_only_icon = "❌"

        table.add_row(
            str(p["id"]),
            p["name"],
            read_only_icon,
        )

    return RichConsoleInterface(retcode=0, data=table)


def get_project(client: AuthClient, *, project_id: int) -> RichJsonInterface:
    res = client.get(f"{settings.BASE_URL}/project/{project_id}/")
    project = check_response(res, expected_status_code=200)
    return RichJsonInterface(retcode=0, data=project)


def delete_project(client: AuthClient, *, project_id: int) -> PrintInterface:

    res = client.delete(f"{settings.BASE_URL}/project/{project_id}/")
    check_response(res, expected_status_code=204)
    return PrintInterface(retcode=0, data="")


def patch_project(
    client: AuthClient,
    *,
    project_id: int,
    new_name: Optional[str] = None,
    make_read_only: bool = False,
    remove_read_only: bool = False,
) -> Union[RichJsonInterface, PrintInterface]:
    project_update = {}
    if new_name:
        project_update["name"] = new_name
    if make_read_only:
        project_update["read_only"] = True
    if remove_read_only:
        project_update["read_only"] = False

    if not project_update:
        return PrintInterface(retcode=1, data="Nothing to update")

    res = client.patch(
        f"{settings.BASE_URL}/project/{project_id}/", json=project_update
    )
    new_project = check_response(res, expected_status_code=200)
    return RichJsonInterface(retcode=0, data=new_project)
