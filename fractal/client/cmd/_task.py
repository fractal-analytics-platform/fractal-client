import json
from typing import Optional

from ..authclient import AuthClient
from ..config import settings
from ..interface import RichJsonInterface
from ..response import check_response
from fractal.common.models import ApplyWorkflow
from fractal.common.models import TaskCreate
from fractal.common.models import TaskRead


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
    name: str,
    resource_type: str,
    input_type: str,
    output_type: str,
    module: str,
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
    return RichJsonInterface(retcode=0, data=new_task.dict())


async def task_edit():
    pass


async def task_add_subtask():
    pass


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
