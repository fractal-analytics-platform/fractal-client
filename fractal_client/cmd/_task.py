import json
import logging
import sys
from typing import Optional

from ..authclient import AuthClient
from ..config import settings
from ..interface import Interface
from ..response import check_response
from ._aux_task_caching import FractalCacheError
from ._aux_task_caching import get_task_id_from_cache
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


def task_collect_custom(
    client: AuthClient,
    *,
    source: str,
    python_interpreter: str,
    manifest: str,
    version: Optional[str] = None,
    package_name: Optional[str] = None,
    package_root: Optional[str] = None,
    batch: bool = False,
) -> Interface:

    try:
        with open(manifest, "r") as f:
            manifest_dict = json.load(f)
    except FileNotFoundError as e:
        raise FileNotFoundError(
            f"Fractal Client cannot find the file {manifest}. "
            "Note that the file must be on the same machine where Fractal "
            f"Client is running.\nOriginal error: {e}."
        )

    task_collect = dict(
        source=source,
        python_interpreter=python_interpreter,
        manifest=manifest_dict,
    )
    if version:
        task_collect["version"] = version
    if package_name:
        task_collect["package_name"] = package_name
    if package_root:
        task_collect["package_root"] = package_root

    res = client.post(
        f"{settings.BASE_URL}/task/collect/custom/", json=task_collect
    )

    task_list = check_response(
        res, expected_status_code=201, redact_long_payload=True
    )

    if batch:
        task_ids = [str(task["id"]) for task in task_list]
        return Interface(retcode=0, data=" ".join(task_ids))
    else:
        return Interface(retcode=0, data=task_list)


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


def post_task(
    client: AuthClient,
    *,
    name: str,
    source: str,
    batch: bool = False,
    command_non_parallel: Optional[str] = None,
    command_parallel: Optional[str] = None,
    version: Optional[str] = None,
    meta_non_parallel: Optional[str] = None,
    meta_parallel: Optional[str] = None,
    args_schema_non_parallel: Optional[str] = None,
    args_schema_parallel: Optional[str] = None,
    args_schema_version: Optional[str] = None,
) -> Interface:
    task = dict(
        name=name,
        source=source,
    )
    if command_non_parallel:
        task["command_non_parallel"] = command_non_parallel
    if command_parallel:
        task["command_parallel"] = command_parallel
    if version:
        task["version"] = version
    if meta_non_parallel:
        with open(meta_non_parallel, "r") as f:
            task["meta_non_parallel"] = json.load(f)
    if meta_parallel:
        with open(meta_parallel, "r") as f:
            task["meta_parallel"] = json.load(f)
    if args_schema_parallel:
        with open(args_schema_parallel, "r") as f:
            task["args_schema_parallel"] = json.load(f)
    if args_schema_non_parallel:
        with open(args_schema_non_parallel, "r") as f:
            task["args_schema_non_parallel"] = json.load(f)
    if args_schema_version:
        task["args_schema_version"] = args_schema_version

    res = client.post(f"{settings.BASE_URL}/task/", json=task)
    new_task = check_response(res, expected_status_code=201)

    if batch:
        return Interface(retcode=0, data=str(new_task["id"]))
    else:
        return Interface(retcode=0, data=new_task)


def patch_task(
    client: AuthClient,
    *,
    id: Optional[int] = None,
    name: Optional[str] = None,
    version: Optional[str] = None,
    new_name: Optional[str] = None,
    new_version: Optional[str] = None,
    command_non_parallel: Optional[str] = None,
    command_parallel: Optional[str] = None,
    input_types: Optional[str] = None,
    output_types: Optional[str] = None,
) -> Interface:

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
    if new_version:
        task_update["version"] = new_version
    if command_non_parallel:
        task_update["command_non_parallel"] = command_non_parallel
    if command_parallel:
        task_update["command_parallel"] = command_parallel
    if input_types:
        with open(input_types, "r") as f:
            task_update["input_types"] = json.load(f)
    if output_types:
        with open(output_types, "r") as f:
            task_update["output_types"] = json.load(f)

    res = client.patch(f"{settings.BASE_URL}/task/{id}/", json=task_update)
    new_task = check_response(res, expected_status_code=200)
    return Interface(retcode=0, data=new_task)


def delete_task(
    client: AuthClient,
    *,
    id: Optional[int] = None,
    name: Optional[str] = None,
    version: Optional[str] = None,
) -> Interface:

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
    return Interface(retcode=0, data="")
