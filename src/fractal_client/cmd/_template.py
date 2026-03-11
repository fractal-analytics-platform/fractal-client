from ..authclient import AuthClient
from ..interface import Interface
from ..response import check_response


def template_show(client: AuthClient, *, template_id: int):
    res = client.get(f"api/v2/workflow_template/{template_id}/")
    data = check_response(res, expected_status_code=200)
    return Interface(retcode=0, data=data)
