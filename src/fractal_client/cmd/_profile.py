import json
from pathlib import Path

from ..authclient import AuthClient
from ..interface import Interface
from ..response import check_response


def post_profile(
    client: AuthClient,
    *,
    resource_id: int,
    json_file: str,
    batch: bool = False,
) -> Interface:
    # Prepare a ResourceCreate request body
    with Path(json_file).open("r") as f:
        payload = json.load(f)

    # Send API request
    res = client.post(
        f"admin/v2/resource/{resource_id}/profile/", json=payload
    )
    profile = check_response(res, expected_status_code=201)
    if batch:
        return Interface(retcode=0, data=profile["id"])
    else:
        return Interface(retcode=0, data=profile)
