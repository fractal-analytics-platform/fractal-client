import json
import logging
from pathlib import Path
from typing import Optional

from ..authclient import AuthClient
from ..common.schemas import ApplyWorkflowCreate
from ..common.schemas import WorkflowCreate
from ..common.schemas import WorkflowRead
from ..common.schemas import WorkflowTaskCreate
from ..common.schemas import WorkflowTaskRead
from ..common.schemas import WorkflowTaskUpdate
from ..common.schemas import WorkflowUpdate
from ..config import settings
from ..interface import BaseInterface
from ..interface import PrintInterface
from ..interface import RichJsonInterface
from ..response import check_response
from .utils import get_cached_task_by_name


async def workflow_new(
    client: AuthClient,
    name: str,
    project_id: int,
    batch: bool = False,
    **kwargs,
) -> BaseInterface:
    workflow = WorkflowCreate(
        name=name,
        project_id=project_id,
    )
    logging.info(workflow)
    res = await client.post(
        f"{settings.BASE_URL}/workflow/",
        json=workflow.dict(),
    )
    workflow = check_response(
        res, expected_status_code=201, coerce=WorkflowRead
    )
    if batch:
        return PrintInterface(retcode=0, data=workflow.id)
    else:
        return RichJsonInterface(retcode=0, data=workflow.dict())


async def workflow_list(
    client: AuthClient,
    project_id: int,
    batch: bool = False,
    **kwargs,
) -> RichJsonInterface:

    res = await client.get(
        f"{settings.BASE_URL}/project/{project_id}/workflows/"
    )
    workflow_list = check_response(res, expected_status_code=200)
    return RichJsonInterface(retcode=0, data=workflow_list)


async def workflow_delete(
    client: AuthClient,
    *,
    id: int,
    **kwargs,
) -> BaseInterface:
    res = await client.delete(f"{settings.BASE_URL}/workflow/{id}")
    check_response(res, expected_status_code=204)
    return PrintInterface(retcode=0, data="")


async def workflow_show(
    client: AuthClient,
    *,
    id: int,
    **kwargs,
) -> RichJsonInterface:
    res = await client.get(f"{settings.BASE_URL}/workflow/{id}")
    workflow = check_response(res, expected_status_code=200)
    return RichJsonInterface(retcode=0, data=workflow)


async def workflow_add_task(
    client: AuthClient,
    *,
    id: int,
    task_id_or_name: str,
    order: int = None,
    args_file: Optional[str] = None,
    **kwargs,
) -> RichJsonInterface:

    try:
        task_id = int(task_id_or_name)
    except ValueError:
        task_id = await get_cached_task_by_name(task_id_or_name, client)

    workflow_task = WorkflowTaskCreate(task_id=task_id, order=order)
    if args_file:
        with Path(args_file).open("r") as f:
            args = json.load(f)
            workflow_task.args = args

    res = await client.post(
        f"{settings.BASE_URL}/workflow/{id}/add-task/",
        json=workflow_task.dict(),
    )
    workflow_task = check_response(
        res, expected_status_code=201, coerce=WorkflowRead
    )

    return RichJsonInterface(retcode=0, data=workflow_task.dict())


async def workflow_edit_task(
    client: AuthClient,
    *,
    id: int,
    workflow_task_id: int,
    json_file: str,
    **kwargs,
) -> RichJsonInterface:

    with Path(json_file).open("r") as f:
        payload = json.load(f)
    payload_update = WorkflowTaskUpdate(**payload)
    res = await client.patch(
        f"{settings.BASE_URL}/workflow/{id}/edit-task/{workflow_task_id}",
        json=payload_update.dict(exclude_unset=True),
    )
    workflow_task = check_response(
        res, expected_status_code=200, coerce=WorkflowTaskRead
    )

    return RichJsonInterface(retcode=0, data=workflow_task.dict())


async def workflow_remove_task(
    client: AuthClient,
    *,
    id: int,
    workflow_task_id: int,
    **kwargs,
) -> BaseInterface:

    res = await client.delete(
        f"{settings.BASE_URL}/workflow/{id}/rm-task/{workflow_task_id}"
    )
    check_response(res, expected_status_code=204)
    return PrintInterface(retcode=0, data="")


async def workflow_edit(
    client: AuthClient,
    *,
    id: str,
    **workflow_update_dict,
) -> BaseInterface:
    workflow_update = WorkflowUpdate(**workflow_update_dict)
    payload = workflow_update.dict(exclude_unset=True)
    if not payload:
        return PrintInterface(retcode=1, data="Nothing to update")

    res = await client.patch(
        f"{settings.BASE_URL}/workflow/{id}", json=payload
    )
    new_workflow = check_response(
        res, expected_status_code=200, coerce=WorkflowRead
    )
    return RichJsonInterface(retcode=0, data=new_workflow.dict())


async def workflow_apply(
    client: AuthClient,
    *,
    workflow_id: int,
    input_dataset_id: int,
    output_dataset_id: int,
    # output_dataset_id: Optional[int] = None,
    overwrite_input: bool = False,
    project_id: Optional[int] = None,
    worker_init: Optional[str] = None,
    **kwargs,
) -> BaseInterface:
    apply_wf_create = ApplyWorkflowCreate(
        workflow_id=workflow_id,
        input_dataset_id=input_dataset_id,
        output_dataset_id=output_dataset_id,
        overwrite_input=overwrite_input,
        project_id=project_id,
        worker_init=worker_init,
    )

    res = await client.post(
        f"{settings.BASE_URL}/project/apply/", json=apply_wf_create.dict()
    )
    # TODO check output

    return RichJsonInterface(retcode=0, data=res.json())
