import json
import os
from pathlib import Path
from typing import List
from typing import Optional
from typing import Union

from ..authclient import AuthClient
from ..config import settings
from ..interface import BaseInterface
from ..interface import PrintInterface
from ..interface import RichJsonInterface
from ..response import check_response
from fractal.common.models import ApplyWorkflow
from fractal.common.models import SubtaskCreate
from fractal.common.models import TaskCreate
from fractal.common.models import TaskRead
from fractal.common.models import TaskUpdate


async def get_cached_task_by_name(name: str, client: AuthClient) -> int:
    # Set paths
    cache_dir = str(Path(f"{settings.FRACTAL_CACHE_PATH}").expanduser())
    cache_file = f"{cache_dir}/tasks"

    # If cache is missing, create it
    cache_up_to_date = False
    if not os.path.isfile(cache_file):
        await refresh_task_cache(client)
        cache_up_to_date = True

    # Read cache
    with open(cache_file, "r") as f:
        task_cache = json.load(f)

    # Look for name in cache
    # Case 1: name exists in cache
    if name in task_cache.keys():
        return task_cache[name]
    # Case 2: name is missing, and cache was just updated
    elif cache_up_to_date:
        raise KeyError(f"Task {name} not in {cache_file}\n")
    # Case 3: name is missing but cache may be out of date
    else:
        await refresh_task_cache(client)
        with open(cache_file, "r") as f:
            task_cache = json.load(f)
        try:
            return task_cache[name]
        except KeyError as e:
            raise KeyError(f"Task {name} not in {cache_file}\n", str(e))


async def refresh_task_cache(client: AuthClient, **kwargs) -> List[dict]:

    # Get task_list
    res = await client.get(f"{settings.BASE_URL}/task/")
    task_list = check_response(res, expected_status_code=200)

    # Refresh cache of (name,id) pairs
    task_cache = {task["name"]: task["id"] for task in task_list}

    # Create cache folder, if needed
    cache_dir = str(Path(f"{settings.FRACTAL_CACHE_PATH}").expanduser())
    if not os.path.isdir(cache_dir):
        os.makedirs(cache_dir)

    # Write task cache
    with open(f"{cache_dir}/tasks", "w") as f:
        json.dump(task_cache, f)

    return task_list


async def task_list(client: AuthClient, **kwargs) -> RichJsonInterface:
    task_list = await refresh_task_cache(client=client, **kwargs)
    return RichJsonInterface(retcode=0, data=task_list)


async def task_new(
    client: AuthClient,
    *,
    batch: bool = False,
    name: str,
    resource_type: str,
    input_type: str,
    output_type: str,
    module: Optional[str] = None,
    default_args: Optional[str] = None,
    subtask_list: Optional[str] = None,
    **kwargs,
) -> RichJsonInterface:

    if not default_args:
        default_args_dict = {}
    else:
        default_args_dict = json.loads(default_args)
    if not subtask_list:
        subtask_list_list = []
    else:
        subtask_list_list = json.loads(subtask_list)

    resource_type = resource_type.replace("_", " ")
    task = TaskCreate(
        name=name,
        resource_type=resource_type,
        input_type=input_type,
        output_type=output_type,
        default_args=default_args_dict,
        module=module,
        subtask_list=subtask_list_list,
    )

    res = await client.post(
        f"{settings.BASE_URL}/task/",
        json=task.dict(),
    )
    new_task = check_response(res, expected_status_code=201, coerce=TaskRead)
    if batch:
        return PrintInterface(retcode=0, data=new_task.id)
    else:
        return RichJsonInterface(retcode=0, data=new_task.dict())


async def task_edit(
    client: AuthClient,
    *,
    task_id_or_name: Union[int, str],
    **task_update_dict,
) -> BaseInterface:
    task_update = TaskUpdate(**task_update_dict)
    payload = task_update.dict(exclude_unset=True)
    if not payload:
        return PrintInterface(retcode=1, data="Nothing to update")

    # FIXME: get ID or name
    task_id = await get_cached_task_by_name(name=task_name, client=client)

    res = await client.patch(
        f"{settings.BASE_URL}/task/{task_id}", json=payload
    )
    new_task = check_response(res, expected_status_code=200, coerce=TaskRead)
    return RichJsonInterface(retcode=0, data=new_task.dict())


async def task_add_subtask(
    client: AuthClient,
    *,
    batch: bool = False,
    parent_task_id_or_name: Union[int, str],
    subtask_name: str,
    args_file: Optional[str] = None,
    order: Optional[int] = None,
    **kwargs,
):

    if args_file:
        with open(args_file, "r") as fin:
            args = json.load(fin)
    else:
        args = {}

    # FIXME: get ID or name
    parent_task_id = await get_cached_task_by_name(
        name=parent_task_name, client=client
    )
    subtask_id = await get_cached_task_by_name(
        name=subtask_name, client=client
    )

    subtask_create = SubtaskCreate(
        parent_task_id=parent_task_id,
        subtask_id=subtask_id,
        order=order,
        args=args,
    )

    res = await client.post(
        f"{settings.BASE_URL}/task/{parent_task_id}/subtask/",
        json=subtask_create.dict(exclude_unset=True),
    )
    new_subtask = check_response(
        res, expected_status_code=201, coerce=TaskRead
    )
    if batch:
        return PrintInterface(retcode=0, data=new_subtask.id)
    else:
        return RichJsonInterface(retcode=0, data=new_subtask.dict())


async def task_apply(
    client: AuthClient,
    *,
    project_id: int,
    input_dataset_id: int,
    output_dataset_id: int,
    workflow_id_or_name: Union[int, str],
    overwrite_input: bool,
    **kwargs,
) -> RichJsonInterface:

    # FIXME: get ID or name
    workflow_id = await get_cached_task_by_name(
        name=workflow_name, client=client
    )

    workflow = ApplyWorkflow(
        project_id=project_id,
        input_dataset_id=input_dataset_id,
        output_dataset_id=output_dataset_id,
        workflow_id=workflow_id,
        overwrite_input=overwrite_input,
    )

    res = await client.post(
        f"{settings.BASE_URL}/project/apply/",
        json=workflow.dict(),
    )
    wf_submitted = check_response(res, expected_status_code=202)
    return RichJsonInterface(retcode=0, data=wf_submitted)
