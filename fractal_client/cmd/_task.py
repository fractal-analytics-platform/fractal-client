import json
import logging
import sys

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
    package_version: str | None = None,
    python_version: str | None = None,
    package_extras: str | None = None,
    pinned_dependency: list[str] | None = None,
    private: bool = False,
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

    is_private = "?private=true" if private else ""

    res = client.post(
        f"{settings.BASE_URL}/task/collect/pip/{is_private}", json=task_collect
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
    label: str,
    python_interpreter: str,
    manifest: str,
    version: str | None = None,
    package_name: str | None = None,
    package_root: str | None = None,
    private: bool = False,
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
        label=label,
        python_interpreter=python_interpreter,
        manifest=manifest_dict,
    )
    if version:
        task_collect["version"] = version
    if package_name:
        task_collect["package_name"] = package_name
    if package_root:
        task_collect["package_root"] = package_root
    is_private = "?private=true" if private else ""

    res = client.post(
        f"{settings.BASE_URL}/task/collect/custom/{is_private}",
        json=task_collect,
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

    res = client.get(f"{settings.BASE_URL}/task/collect/{state_id}/")
    state = check_response(res, expected_status_code=200)

    # Remove key-value pairs with None value
    state["data"] = {key: val for (key, val) in state["data"].items() if val}
    if (include_logs is False) and ("log" in state["data"]):
        state["data"]["log"] = None

    return Interface(retcode=0, data=state)


def post_task(
    client: AuthClient,
    *,
    name: str,
    batch: bool = False,
    command_non_parallel: str | None = None,
    command_parallel: str | None = None,
    version: str | None = None,
    meta_non_parallel: str | None = None,
    meta_parallel: str | None = None,
    args_schema_non_parallel: str | None = None,
    args_schema_parallel: str | None = None,
    args_schema_version: str | None = None,
    private: bool = False,
) -> Interface:
    task = dict(name=name)
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
    is_private = "?private=true" if private else ""

    res = client.post(f"{settings.BASE_URL}/task/{is_private}", json=task)
    new_task = check_response(res, expected_status_code=201)

    if batch:
        return Interface(retcode=0, data=str(new_task["id"]))
    else:
        return Interface(retcode=0, data=new_task)


def patch_task(
    client: AuthClient,
    *,
    id: int | None = None,
    name: str | None = None,
    version: str | None = None,
    command_non_parallel: str | None = None,
    command_parallel: str | None = None,
    input_types: str | None = None,
    output_types: str | None = None,
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
