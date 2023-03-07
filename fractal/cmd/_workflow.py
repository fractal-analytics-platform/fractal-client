import json
import logging
from pathlib import Path
from typing import Optional

from ..authclient import AuthClient
from ..common.schemas import ApplyWorkflowCreate
from ..common.schemas import ApplyWorkflowRead
from ..common.schemas import WorkflowCreate
from ..common.schemas import WorkflowExport
from ..common.schemas import WorkflowImport
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
    workflow_id: int,
    **kwargs,
) -> BaseInterface:
    res = await client.delete(f"{settings.BASE_URL}/workflow/{workflow_id}")
    check_response(res, expected_status_code=204)
    return PrintInterface(retcode=0, data="")


async def workflow_show(
    client: AuthClient,
    *,
    workflow_id: int,
    **kwargs,
) -> RichJsonInterface:
    res = await client.get(f"{settings.BASE_URL}/workflow/{workflow_id}")
    workflow = check_response(
        res, expected_status_code=200, coerce=WorkflowRead
    )
    return RichJsonInterface(retcode=0, data=workflow.dict())


async def workflow_add_task(
    client: AuthClient,
    *,
    workflow_id: int,
    batch: bool = False,
    task_id_or_name: str,
    order: Optional[int] = None,
    args_file: Optional[str] = None,
    meta_file: Optional[str] = None,
    **kwargs,
) -> RichJsonInterface:

    try:
        task_id = int(task_id_or_name)
    except ValueError:
        task_id = await get_cached_task_by_name(task_id_or_name, client)

    if order is None:
        workflow_task = WorkflowTaskCreate(task_id=task_id)
    else:
        workflow_task = WorkflowTaskCreate(task_id=task_id, order=order)
    if args_file:
        with Path(args_file).open("r") as f:
            args = json.load(f)
            workflow_task.args = args
    if meta_file:
        with Path(meta_file).open("r") as f:
            meta = json.load(f)
            workflow_task.meta = meta

    res = await client.post(
        f"{settings.BASE_URL}/workflow/{workflow_id}/add-task/",
        json=workflow_task.dict(exclude_unset=True),
    )
    workflow_task = check_response(
        res, expected_status_code=201, coerce=WorkflowTaskRead
    )

    if batch:
        return PrintInterface(retcode=0, data=str(workflow_task.id))
    else:
        return RichJsonInterface(retcode=0, data=workflow_task.dict())


async def workflow_edit_task(
    client: AuthClient,
    *,
    workflow_id: int,
    workflow_task_id: int,
    args_file: Optional[str] = None,
    meta_file: Optional[str] = None,
    **kwargs,
) -> RichJsonInterface:

    # Check that at least one of args_file or meta_file was given (note: it
    # would be reasonable to check it in the parser, but we are not aware of a
    # method within argparse).
    if not (args_file or meta_file):
        raise ValueError(
            "At least one of {args_file,meta_file} arguments is required"
        )

    # Combine args/meta payload
    payload = {}
    if args_file:
        with Path(args_file).open("r") as f:
            args = json.load(f)
            payload["args"] = args
    if meta_file:
        with Path(meta_file).open("r") as f:
            meta = json.load(f)
            payload["meta"] = meta

    payload_update = WorkflowTaskUpdate(**payload)
    res = await client.patch(
        f"{settings.BASE_URL}/"
        f"workflow/{workflow_id}/edit-task/{workflow_task_id}",
        json=payload_update.dict(exclude_unset=True),
    )
    workflow_task = check_response(
        res, expected_status_code=200, coerce=WorkflowTaskRead
    )

    return RichJsonInterface(retcode=0, data=workflow_task.dict())


async def workflow_remove_task(
    client: AuthClient,
    *,
    workflow_id: int,
    workflow_task_id: int,
    **kwargs,
) -> BaseInterface:

    res = await client.delete(
        f"{settings.BASE_URL}/"
        f"workflow/{workflow_id}/rm-task/{workflow_task_id}"
    )
    check_response(res, expected_status_code=204)
    return PrintInterface(retcode=0, data="")


async def workflow_edit(
    client: AuthClient,
    *,
    workflow_id: int,
    **workflow_update_dict,
) -> BaseInterface:
    workflow_update = WorkflowUpdate(**workflow_update_dict)
    payload = workflow_update.dict(exclude_unset=True)
    if not payload:
        return PrintInterface(retcode=1, data="Nothing to update")

    res = await client.patch(
        f"{settings.BASE_URL}/workflow/{workflow_id}", json=payload
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
    apply_wf_create_dict = dict(
        workflow_id=workflow_id,
        input_dataset_id=input_dataset_id,
        output_dataset_id=output_dataset_id,
        overwrite_input=overwrite_input,
    )
    # Prepare ApplyWorkflowCreate object, without None attributes
    if project_id:
        apply_wf_create_dict["project_id"] = project_id
    if worker_init:
        apply_wf_create_dict["worker_init"] = worker_init
    apply_wf_create = ApplyWorkflowCreate(**apply_wf_create_dict)

    res = await client.post(
        f"{settings.BASE_URL}/project/apply/",
        json=apply_wf_create.dict(exclude_unset=True),
    )
    apply_wf_read = check_response(
        res, expected_status_code=202, coerce=ApplyWorkflowRead
    )
    return RichJsonInterface(retcode=0, data=apply_wf_read.sanitised_dict())


async def workflow_import(
    client: AuthClient,
    *,
    project_id: int,
    json_file: str,
    batch: bool = False,
    **kwargs,
) -> BaseInterface:
    with Path(json_file).open("r") as f:
        workflow = json.load(f)
    workflow = WorkflowImport(**workflow)
    res = await client.post(
        f"{settings.BASE_URL}/project/{project_id}/import-workflow/",
        json=workflow.dict(exclude_unset=True),
    )
    wf_read = check_response(
        res, expected_status_code=201, coerce=WorkflowRead
    )
    if batch:
        datastr = f"{wf_read.id}"
        for wftask in wf_read.task_list:
            datastr += f" {wftask.id}"
        return PrintInterface(retcode=0, data=datastr)
    else:
        return RichJsonInterface(retcode=0, data=wf_read.dict())


async def workflow_export(
    client: AuthClient,
    *,
    workflow_id: int,
    json_file: str,
    **kwargs,
) -> BaseInterface:
    res = await client.get(
        f"{settings.BASE_URL}/workflow/{workflow_id}/export/"
    )
    workflow = check_response(
        res, expected_status_code=200, coerce=WorkflowExport
    )
    with Path(json_file).open("w") as f:
        json.dump(workflow.dict(), f, indent=2)
    return PrintInterface(
        retcode=0, data=f"Workflow {workflow_id} exported at {json_file}"
    )
