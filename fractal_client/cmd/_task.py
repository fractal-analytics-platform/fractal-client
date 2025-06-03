import json
import logging
import sys

from ..authclient import AuthClient
from ..interface import Interface
from ..response import check_response
from ._aux_task_caching import FractalCacheError
from ._aux_task_caching import get_task_id_from_cache
from ._aux_task_caching import refresh_task_cache


def get_task_list(client: AuthClient) -> Interface:
    task_list = refresh_task_cache(client=client)
    return Interface(retcode=0, data=task_list)


def post_task(
    client: AuthClient,
    *,
    name: str,
    task_type: str | None = None,
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
    if task_type:
        task["type"] = task_type
    if command_non_parallel:
        task["command_non_parallel"] = command_non_parallel
    if command_parallel:
        task["command_parallel"] = command_parallel
    if version:
        task["version"] = version
    if meta_non_parallel:
        with open(meta_non_parallel) as f:
            task["meta_non_parallel"] = json.load(f)
    if meta_parallel:
        with open(meta_parallel) as f:
            task["meta_parallel"] = json.load(f)
    if args_schema_parallel:
        with open(args_schema_parallel) as f:
            task["args_schema_parallel"] = json.load(f)
    if args_schema_non_parallel:
        with open(args_schema_non_parallel) as f:
            task["args_schema_non_parallel"] = json.load(f)
    if args_schema_version:
        task["args_schema_version"] = args_schema_version
    is_private = "?private=true" if private else ""

    res = client.post(f"api/v2/task/{is_private}", json=task)
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
        with open(input_types) as f:
            task_update["input_types"] = json.load(f)
    if output_types:
        with open(output_types) as f:
            task_update["output_types"] = json.load(f)

    res = client.patch(f"api/v2/task/{id}/", json=task_update)
    new_task = check_response(res, expected_status_code=200)
    return Interface(retcode=0, data=new_task)
