import json
from typing import Optional

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


async def task_list(
    client: AuthClient,
    **kwargs,
) -> RichJsonInterface:
    res = await client.get(
        f"{settings.BASE_URL}/task/",
    )
    task_list = check_response(res, expected_status_code=200)
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
    task_id: int,
    **task_update_dict,
) -> BaseInterface:
    task_update = TaskUpdate(**task_update_dict)
    payload = task_update.dict(exclude_unset=True)
    if not payload:
        return PrintInterface(retcode=1, data="Nothing to update")

    res = await client.patch(
        f"{settings.BASE_URL}/task/{task_id}", json=payload
    )
    new_task = check_response(res, expected_status_code=200, coerce=TaskRead)
    return RichJsonInterface(retcode=0, data=new_task.dict())


async def task_add_subtask(
    client: AuthClient,
    *,
    batch: bool = False,
    parent_task_id: int,
    subtask_id: int,
    args_file: Optional[str] = None,
    order: Optional[int] = None,
    **kwargs,
):

    if args_file:
        with open(args_file, "r") as fin:
            args = json.load(fin)
    else:
        args = {}

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
    workflow_id: int,
    overwrite_input: bool,
    **kwargs,
) -> RichJsonInterface:

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
