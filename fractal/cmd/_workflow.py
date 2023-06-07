import json
import logging
import sys
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
from ._aux_task_caching import FractalCacheError
from ._aux_task_caching import get_task_id_from_cache


async def post_workflow(
    client: AuthClient, *, name: str, project_id: int, batch: bool = False
) -> BaseInterface:
    workflow = WorkflowCreate(
        name=name,
        project_id=project_id,
    )
    res = await client.post(
        f"{settings.BASE_URL}/project/{project_id}/workflow/",
        json=workflow.dict(),
    )
    workflow = check_response(
        res, expected_status_code=201, coerce=WorkflowRead
    )
    if batch:
        return PrintInterface(retcode=0, data=workflow.id)
    else:
        return RichJsonInterface(retcode=0, data=workflow.dict())


async def get_workflow_list(
    client: AuthClient, *, project_id: int, batch: bool = False
) -> RichJsonInterface:

    res = await client.get(
        f"{settings.BASE_URL}/project/{project_id}/workflow/"
    )
    workflow_list = check_response(res, expected_status_code=200)
    return RichJsonInterface(retcode=0, data=workflow_list)


async def delete_workflow(
    client: AuthClient, *, project_id: int, workflow_id: int
) -> BaseInterface:
    res = await client.delete(
        f"{settings.BASE_URL}/project/{project_id}/workflow/{workflow_id}"
    )
    check_response(res, expected_status_code=204)
    return PrintInterface(retcode=0, data="")


async def get_workflow(
    client: AuthClient, *, project_id: int, workflow_id: int
) -> RichJsonInterface:
    res = await client.get(
        f"{settings.BASE_URL}/project/{project_id}/workflow/{workflow_id}"
    )
    workflow = check_response(
        res, expected_status_code=200, coerce=WorkflowRead
    )
    return RichJsonInterface(retcode=0, data=workflow.dict())


async def post_workflowtask(
    client: AuthClient,
    *,
    project_id: int,
    workflow_id: int,
    task_id: Optional[int] = None,
    task_name: Optional[str] = None,
    batch: bool = False,
    version: Optional[str] = None,
    order: Optional[int] = None,
    args_file: Optional[str] = None,
    meta_file: Optional[str] = None,
) -> RichJsonInterface:

    if task_id and version:
        logging.warning(
            "Too many arguments: cannot provide both "
            "`task_id` and `task_version`."
        )
        sys.exit(1)
    elif task_name:
        try:
            task_id = await get_task_id_from_cache(
                client=client, task_name=task_name, version=version
            )
        except FractalCacheError as e:
            print(e)
            sys.exit(1)

    if order is None:
        workflow_task = WorkflowTaskCreate()
    else:
        workflow_task = WorkflowTaskCreate(order=order)
    if args_file:
        with Path(args_file).open("r") as f:
            args = json.load(f)
            workflow_task.args = args
    if meta_file:
        with Path(meta_file).open("r") as f:
            meta = json.load(f)
            workflow_task.meta = meta

    res = await client.post(
        (
            f"{settings.BASE_URL}/project/{project_id}/"
            f"workflow/{workflow_id}/wftask/"
            f"?{task_id=}"
        ),
        json=workflow_task.dict(exclude_unset=True),
    )
    workflow_task = check_response(
        res, expected_status_code=201, coerce=WorkflowTaskRead
    )

    if batch:
        return PrintInterface(retcode=0, data=str(workflow_task.id))
    else:
        return RichJsonInterface(retcode=0, data=workflow_task.dict())


async def patch_workflowtask(
    client: AuthClient,
    *,
    project_id: int,
    workflow_id: int,
    workflow_task_id: int,
    args_file: Optional[str] = None,
    meta_file: Optional[str] = None,
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
        (
            f"{settings.BASE_URL}/project/{project_id}/"
            f"workflow/{workflow_id}/wftask/{workflow_task_id}"
        ),
        json=payload_update.dict(exclude_unset=True),
    )
    workflow_task = check_response(
        res, expected_status_code=200, coerce=WorkflowTaskRead
    )

    return RichJsonInterface(retcode=0, data=workflow_task.dict())


async def delete_workflowtask(
    client: AuthClient,
    *,
    project_id: int,
    workflow_id: int,
    workflow_task_id: int,
) -> BaseInterface:

    res = await client.delete(
        f"{settings.BASE_URL}/project/{project_id}/"
        f"workflow/{workflow_id}/wftask/{workflow_task_id}"
    )
    check_response(res, expected_status_code=204)
    return PrintInterface(retcode=0, data="")


async def patch_workflow(
    client: AuthClient,
    *,
    project_id: int,
    workflow_id: int,
    new_name: str,
) -> BaseInterface:

    workflow_update = WorkflowUpdate(name=new_name)
    payload = workflow_update.dict(exclude_unset=True)

    res = await client.patch(
        f"{settings.BASE_URL}/project/{project_id}/workflow/{workflow_id}",
        json=payload,
    )
    new_workflow = check_response(
        res, expected_status_code=200, coerce=WorkflowRead
    )
    return RichJsonInterface(retcode=0, data=new_workflow.dict())


async def workflow_apply(
    client: AuthClient,
    *,
    project_id: int,
    workflow_id: int,
    input_dataset_id: int,
    output_dataset_id: int,
    worker_init: Optional[str] = None,
) -> BaseInterface:
    apply_wf_create_dict = dict(
        workflow_id=workflow_id,
        input_dataset_id=input_dataset_id,
        output_dataset_id=output_dataset_id,
    )
    # Prepare ApplyWorkflowCreate object, without None attributes
    if worker_init:
        apply_wf_create_dict["worker_init"] = worker_init
    apply_wf_create = ApplyWorkflowCreate(**apply_wf_create_dict)

    # Prepare query parameters
    query_parameters = (
        f"input_dataset_id={input_dataset_id}"
        f"&output_dataset_id={output_dataset_id}"
    )

    res = await client.post(
        (
            f"{settings.BASE_URL}/project/{project_id}/workflow/{workflow_id}/"
            f"apply/?{query_parameters}"
        ),
        json=apply_wf_create.dict(exclude_unset=True),
    )
    apply_wf_read = check_response(
        res, expected_status_code=202, coerce=ApplyWorkflowRead
    )
    return RichJsonInterface(retcode=0, data=apply_wf_read.sanitised_dict())


async def workflow_import(
    client: AuthClient, *, project_id: int, json_file: str, batch: bool = False
) -> BaseInterface:
    with Path(json_file).open("r") as f:
        workflow = json.load(f)
    workflow = WorkflowImport(**workflow)
    res = await client.post(
        f"{settings.BASE_URL}/project/{project_id}/workflow/import/",
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
    project_id: int,
    workflow_id: int,
    json_file: str,
) -> BaseInterface:
    res = await client.get(
        (
            f"{settings.BASE_URL}/project/{project_id}/"
            f"workflow/{workflow_id}/export/"
        ),
    )
    workflow = check_response(
        res, expected_status_code=200, coerce=WorkflowExport
    )
    with Path(json_file).open("w") as f:
        json.dump(workflow.dict(), f, indent=2)
    return PrintInterface(
        retcode=0, data=f"Workflow {workflow_id} exported at {json_file}"
    )
