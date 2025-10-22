import json
from pathlib import Path
from ..authclient import AuthClient
from ..interface import Interface
from ..response import check_response


def post_resource(
    client: AuthClient,
    *,
    json_file: str,
    batch: bool = False,
) -> Interface:
    # Prepare a ResourceCreate request body
    with Path(json_file).open("r") as f:
        payload = json.load(f)

    # Send API request
    res = client.post("admin/v2/resource/", json=payload)
    resource = check_response(res, expected_status_code=201)
    if batch:
        return Interface(retcode=0, data=resource["id"])
    else:
        return Interface(retcode=0, data=resource)
