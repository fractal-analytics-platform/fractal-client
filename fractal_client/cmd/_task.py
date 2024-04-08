import logging
import sys
from typing import Optional

from ..authclient import AuthClient
from ..config import settings
from ..interface import Interface
from ..response import check_response
from ._aux_task_caching import refresh_task_cache


def get_task_list(client: AuthClient) -> Interface:
    task_list = refresh_task_cache(client=client)
    return Interface(retcode=0, data=task_list)


def task_collect_pip(
    client: AuthClient,
    *,
    package: str,
    package_version: Optional[str] = None,
    python_version: Optional[str] = None,
    package_extras: Optional[str] = None,
    pinned_dependency: Optional[list[str]] = None,
    batch: bool = False,
) -> Interface:

    # Construct TaskCollectPip object
    task_collect = dict(package=package)
    if package_version:
        task_collect["package_version"] = package_version
    if python_version:
        task_collect["python_version"] = python_version
    if package_extras:
        task_collect["package_extras"] = package_extras
    if pinned_dependency:
        for pin in pinned_dependency:
            if len(pin.split("=")) != 2:
                logging.error(
                    f"Invalid pin: {pin}.\nPins must be written as "
                    "'--pinned-dependency PACKAGE_NAME=PACKAGE_VERSION'"
                )
                sys.exit(1)
        task_collect["pinned_package_versions"] = {
            _name: _version
            for _name, _version in (p.split("=") for p in pinned_dependency)
        }

    res = client.post(
        f"{settings.BASE_URL}/task/collect/pip/", json=task_collect
    )

    state = check_response(res, expected_status_code=[200, 201])
    if batch:
        output = f"{state['id']} {state['data']['venv_path']}"
        return Interface(retcode=0, data=output)
    else:
        return Interface(retcode=0, data=state)


def task_collection_check(
    client: AuthClient,
    *,
    state_id: int,
    include_logs: bool,
) -> Interface:

    res = client.get(
        f"{settings.BASE_URL}/task/collect/{state_id}/?verbose={include_logs}"
    )
    state = check_response(res, expected_status_code=200)

    # Remove key-value pairs with None value
    state["data"] = {key: val for (key, val) in state["data"].items() if val}

    return Interface(retcode=0, data=state)


def post_task() -> None:
    raise NotImplementedError


def patch_task() -> None:
    raise NotImplementedError


def delete_task() -> None:
    raise NotImplementedError
