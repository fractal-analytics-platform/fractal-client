import json
import logging
import sys
from typing import Optional

from ..authclient import AuthClient
from ..config import settings
from ..interface import BaseInterface
from ..interface import PrintInterface
from ..interface import RichJsonInterface
from ..response import check_response
from ._aux_task_caching import FractalCacheError
from ._aux_task_caching import get_task_id_from_cache
from ._aux_task_caching import refresh_task_cache


def get_task_list(client: AuthClient) -> RichJsonInterface:
    task_list = refresh_task_cache(client=client)
    return RichJsonInterface(retcode=0, data=task_list)


def task_collect_pip(
    client: AuthClient,
    *,
    package: str,
    package_version: Optional[str] = None,
    python_version: Optional[str] = None,
    package_extras: Optional[str] = None,
    pinned_dependency: Optional[list[str]] = None,
    batch: bool = False,
) -> BaseInterface:

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
        return PrintInterface(retcode=0, data=output)
    else:
        return RichJsonInterface(retcode=0, data=state)


def task_collection_check(
    client: AuthClient,
    *,
    state_id: int,
    include_logs: bool,
    do_not_separate_logs: bool = False,
) -> BaseInterface:

    res = client.get(
        f"{settings.BASE_URL}/task/collect/{state_id}/?verbose={include_logs}"
    )
    state = check_response(res, expected_status_code=200)

    # Remove key-value pairs with None value
    state["data"] = {key: val for (key, val) in state["data"].items() if val}

    if (not include_logs) or do_not_separate_logs:
        return RichJsonInterface(retcode=0, data=state)
    else:
        log = state["data"].pop("log")
        extra_lines = ["\nThis is the task-collection log:\n", log]
        return RichJsonInterface(
            retcode=0, data=state, extra_lines=extra_lines
        )


def post_task(
    client: AuthClient,
    *,
    name: str,
    command: str,
    source: str,
    input_type: str = "Any",
    output_type: str = "Any",
    batch: bool = False,
    version: Optional[str] = None,
    meta_file: Optional[str] = None,
    args_schema: Optional[str] = None,
    args_schema_version: Optional[str] = None,
) -> BaseInterface:
    task = dict(
        name=name,
        command=command,
        source=source,
        input_type=input_type,
        output_type=output_type,
    )
    if version:
        task["version"] = version
    if meta_file:
        with open(meta_file, "r") as f:
            task["meta"] = json.load(f)
    if args_schema:
        with open(args_schema, "r") as f:
            task["args_schema"] = json.load(f)
    if args_schema_version:
        task["args_schema_version"] = args_schema_version

    res = client.post(f"{settings.BASE_URL}/task/", json=task)
    new_task = check_response(res, expected_status_code=201)

    if batch:
        return PrintInterface(retcode=0, data=str(new_task["id"]))
    else:
        return RichJsonInterface(retcode=0, data=new_task)


def patch_task(
    client: AuthClient,
    *,
    id: Optional[int] = None,
    name: Optional[str] = None,
    version: Optional[str] = None,
    new_name: Optional[str] = None,
    new_command: Optional[str] = None,
    new_input_type: Optional[str] = None,
    new_output_type: Optional[str] = None,
    new_version: Optional[str] = None,
    meta_file: Optional[str] = None,
    new_args_schema: Optional[str] = None,
    new_args_schema_version: Optional[str] = None,
) -> BaseInterface:

    if id:
        if version:
            logging.error(
                "Too many arguments: cannot provide both `id` and `version`."
            )
            sys.exit(1)
    else:
        try:
            id = get_task_id_from_cache(
                client=client, task_name=name, version=version
            )
        except FractalCacheError as e:
            print(e)
            sys.exit(1)

    task_update = {}
    if new_name:
        task_update["name"] = new_name
    if new_input_type:
        task_update["input_type"] = new_input_type
    if new_output_type:
        task_update["output_type"] = new_output_type
    if new_command:
        task_update["command"] = new_command
    if new_version:
        task_update["version"] = new_version
    if meta_file:
        with open(meta_file, "r") as f:
            task_update["meta"] = json.load(f)
    if new_args_schema:
        with open(new_args_schema, "r") as f:
            task_update["args_schema"] = json.load(f)
    if new_args_schema_version:
        task_update["args_schema_version"] = new_args_schema_version

    if not task_update:
        return PrintInterface(retcode=1, data="Nothing to update")

    res = client.patch(f"{settings.BASE_URL}/task/{id}/", json=task_update)
    new_task = check_response(res, expected_status_code=200)
    return RichJsonInterface(retcode=0, data=new_task)


def delete_task(
    client: AuthClient,
    *,
    id: Optional[int] = None,
    name: Optional[str] = None,
    version: Optional[str] = None,
) -> PrintInterface:

    if id:
        if version:
            logging.error(
                "Too many arguments: cannot provide both `id` and `version`."
            )
            sys.exit(1)
    else:
        try:
            id = get_task_id_from_cache(
                client=client, task_name=name, version=version
            )
        except FractalCacheError as e:
            print(e)
            sys.exit(1)
    res = client.delete(f"{settings.BASE_URL}/task/{id}/")
    check_response(res, expected_status_code=204)
    return PrintInterface(retcode=0, data="")
