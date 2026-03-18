import json
import sys
from pathlib import Path

from ..authclient import AuthClient
from ..interface import Interface
from ..response import check_response


def template_show(client: AuthClient, *, template_id: int):
    res = client.get(f"api/v2/workflow-template/{template_id}/")
    data = check_response(res, expected_status_code=200)
    return Interface(retcode=0, data=data)


def template_new(
    client: AuthClient,
    *,
    workflow_id: int | None = None,
    json_file: str | None = None,
    name: str | None = None,
    version: int | None = None,
    user_group_id: int | None = None,
    batch: bool = False,
):
    if workflow_id is not None:
        if name is None or version is None:
            print(
                "You must provide '--name' and '--version' when you create a "
                "template from a workflow."
            )
            sys.exit(1)
        payload = {
            "name": name,
            "version": version,
        }
        url = f"api/v2/workflow-template/?workflow_id={workflow_id}"
        if user_group_id is not None:
            url += f"&user_group_id={user_group_id}"
    else:  # json_file is not None
        with open(json_file, "r") as f:
            payload = json.load(f)
        if name is not None:
            payload["name"] = name
        if version is not None:
            payload["version"] = version
        url = "api/v2/workflow-template/import/"
        if user_group_id is not None:
            url += f"?user_group_id={user_group_id}"

    res = client.post(url, json=payload)
    template = check_response(res, expected_status_code=201)

    if batch:
        return Interface(retcode=0, data=template["id"])
    else:
        return Interface(retcode=0, data=template)


def template_delete(client: AuthClient, *, template_id: int):
    res = client.delete(f"api/v2/workflow-template/{template_id}/")
    check_response(res, expected_status_code=204)
    return Interface(retcode=0, data="")


def template_export(
    client: AuthClient,
    *,
    template_id: int,
    json_file: str,
) -> Interface:
    res = client.get(
        (f"api/v2/workflow-template/{template_id}/export/"),
    )
    template = check_response(res, expected_status_code=200)

    # mode="x" means "open for exclusive creation, failing if the file
    # already exists" https://docs.python.org/3/library/functions.html#open
    with Path(json_file).open(mode="x") as f:
        json.dump(template, f, indent=2)
    return Interface(
        retcode=0, data=f"Template {template_id} exported at {json_file}."
    )
