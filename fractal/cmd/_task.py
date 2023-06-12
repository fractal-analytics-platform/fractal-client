import json
import logging
import sys
from typing import Optional

from ..authclient import AuthClient
from ..common.schemas import StateRead
from ..common.schemas import TaskCollectPip
from ..common.schemas import TaskCreate
from ..common.schemas import TaskRead
from ..common.schemas import TaskUpdate
from ..config import settings
from ..interface import BaseInterface
from ..interface import PrintInterface
from ..interface import RichJsonInterface
from ..response import check_response
from ._aux_task_caching import FractalCacheError
from ._aux_task_caching import get_task_id_from_cache
from ._aux_task_caching import refresh_task_cache


async def get_task_list(client: AuthClient) -> RichJsonInterface:
    task_list = await refresh_task_cache(client=client)
    return RichJsonInterface(retcode=0, data=task_list)


async def task_collect_pip(
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
    attributes = dict(package=package)
    if package_version:
        attributes["package_version"] = package_version
    if python_version:
        attributes["python_version"] = python_version
    if package_extras:
        attributes["package_extras"] = package_extras
    if pinned_dependency:
        for pin in pinned_dependency:
            if len(pin.split("=")) != 2:
                logging.error(
                    f"Invalid pin: {pin}.\nPins must be written as "
                    "'--pinned-dependency PACKAGE_NAME=PACKAGE_VERSION'"
                )
                sys.exit(1)
        attributes["pinned_package_versions"] = {
            _name: _version
            for _name, _version in (p.split("=") for p in pinned_dependency)
        }
    task_collect = TaskCollectPip(**attributes)

    res = await client.post(
        f"{settings.BASE_URL}/task/collect/pip/",
        json=task_collect.dict(exclude_unset=True),
    )

    state = check_response(
        res, expected_status_code=[200, 201], coerce=StateRead
    )
    if batch:
        output = f"{state.id} {state.data['venv_path']}"
        return PrintInterface(retcode=0, data=output)
    else:
        return RichJsonInterface(retcode=0, data=state.sanitised_dict())


async def task_collection_check(
    client: AuthClient,
    *,
    state_id: int,
    include_logs: bool,
    do_not_separate_logs: bool = False,
) -> BaseInterface:

    res = await client.get(
        f"{settings.BASE_URL}/task/collect/{state_id}?verbose={include_logs}"
    )
    state = check_response(res, expected_status_code=200, coerce=StateRead)

    state_dict = state.sanitised_dict()

    # Remove key-value pairs with None value
    state_dict["data"] = {
        key: val for (key, val) in state_dict["data"].items() if val
    }

    if (not include_logs) or do_not_separate_logs:
        return RichJsonInterface(retcode=0, data=state_dict)
    else:
        log = state_dict["data"].pop("log")
        extra_lines = ["\nThis is the task-collection log:\n", log]
        return RichJsonInterface(
            retcode=0, data=state_dict, extra_lines=extra_lines
        )


async def post_task(
    client: AuthClient,
    *,
    name: str,
    command: str,
    source: str,
    input_type: str = "Any",
    output_type: str = "Any",
    version: Optional[str] = None,
    default_args_file: Optional[str] = None,
    meta_file: Optional[str] = None,
    batch: bool = False,
) -> BaseInterface:
    optionals = {}
    if version:
        optionals["version"] = version
    if default_args_file:
        with open(default_args_file, "r") as f:
            optionals["default_args"] = json.load(f)
    if meta_file:
        with open(meta_file, "r") as f:
            optionals["meta"] = json.load(f)
    payload = TaskCreate(
        name=name,
        command=command,
        source=source,
        input_type=input_type,
        output_type=output_type,
        **optionals,
    ).dict(exclude_unset=True)
    res = await client.post(f"{settings.BASE_URL}/task/", json=payload)
    new_task = check_response(res, expected_status_code=201, coerce=TaskRead)

    if batch:
        return PrintInterface(retcode=0, data=str(new_task.id))
    else:
        return RichJsonInterface(retcode=0, data=new_task.dict())


async def patch_task(
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
) -> BaseInterface:

    if id:
        if version:
            logging.error(
                "Too many arguments: cannot provide both `id` and `version`."
            )
            sys.exit(1)
    else:
        try:
            id = await get_task_id_from_cache(
                client=client, task_name=name, version=version
            )
        except FractalCacheError as e:
            print(e)
            sys.exit(1)

    update = {}
    if new_name:
        update["name"] = new_name
    if new_input_type:
        update["input_type"] = new_input_type
    if new_output_type:
        update["output_type"] = new_output_type
    if new_command:
        update["command"] = new_command
    if new_version:
        update["version"] = new_version
    if meta_file:
        with open(meta_file, "r") as f:
            update["meta"] = json.load(f)

    task_update = TaskUpdate(**update)
    payload = task_update.dict(exclude_unset=True)
    if not payload:
        return PrintInterface(retcode=1, data="Nothing to update")

    res = await client.patch(f"{settings.BASE_URL}/task/{id}", json=payload)
    new_task = check_response(res, expected_status_code=200, coerce=TaskRead)
    return RichJsonInterface(retcode=0, data=new_task.dict())


async def delete_task(
    client: AuthClient,
    *,
    id: Optional[int] = None,
    name: Optional[str] = None,
    version: Optional[str] = None,
) -> PrintInterface:

    raise NotImplementedError("task_delete")
