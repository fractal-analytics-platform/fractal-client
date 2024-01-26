import json
import logging
import sys
from pathlib import Path
from typing import Optional

from ..authclient import AuthClient
from ..config import settings
from ..interface import BaseInterface
from ..interface import PrintInterface
from ..interface import RichJsonInterface
from ..response import check_response
from ._aux_task_caching import FractalCacheError
from ._aux_task_caching import get_task_id_from_cache


def post_workflow(
    client: AuthClient, *, name: str, project_id: int, batch: bool = False
) -> BaseInterface:
    workflow = dict(
        name=name,
        project_id=project_id,
    )
    res = client.post(
        f"{settings.BASE_URL}/project/{project_id}/workflow/",
        json=workflow,
    )
    workflow = check_response(res, expected_status_code=201)
    if batch:
        return PrintInterface(retcode=0, data=workflow["id"])
    else:
        return RichJsonInterface(retcode=0, data=workflow)


def get_workflow_list(
    client: AuthClient, *, project_id: int, batch: bool = False
) -> RichJsonInterface:

    res = client.get(f"{settings.BASE_URL}/project/{project_id}/workflow/")
    workflow_list = check_response(res, expected_status_code=200)
    return RichJsonInterface(retcode=0, data=workflow_list)


def delete_workflow(
    client: AuthClient, *, project_id: int, workflow_id: int
) -> BaseInterface:
    res = client.delete(
        f"{settings.BASE_URL}/project/{project_id}/workflow/{workflow_id}/"
    )
    check_response(res, expected_status_code=204)
    return PrintInterface(retcode=0, data="")


def get_workflow(
    client: AuthClient, *, project_id: int, workflow_id: int
) -> RichJsonInterface:
    res = client.get(
        f"{settings.BASE_URL}/project/{project_id}/workflow/{workflow_id}/"
    )
    workflow = check_response(res, expected_status_code=200)
    return RichJsonInterface(retcode=0, data=workflow)


def post_workflowtask(
    client: AuthClient,
    *,
    project_id: int,
    workflow_id: int,
    task_id: Optional[int] = None,
    task_name: Optional[str] = None,
    task_version: Optional[str] = None,
    batch: bool = False,
    order: Optional[int] = None,
    args_file: Optional[str] = None,
    meta_file: Optional[str] = None,
) -> RichJsonInterface:

    if task_id:
        if task_version:
            logging.error(
                "Too many arguments: cannot provide both "
                "`task_id` and `task_version`."
            )
            sys.exit(1)
    else:
        try:
            task_id = get_task_id_from_cache(
                client=client, task_name=task_name, version=task_version
            )
        except FractalCacheError as e:
            print(e)
            sys.exit(1)

    if order is None:
        workflow_task = dict()
    else:
        workflow_task = dict(order=order)
    if args_file:
        with Path(args_file).open("r") as f:
            args = json.load(f)
            workflow_task["args"] = args
    if meta_file:
        with Path(meta_file).open("r") as f:
            meta = json.load(f)
            workflow_task["meta"] = meta

    res = client.post(
        (
            f"{settings.BASE_URL}/project/{project_id}/"
            f"workflow/{workflow_id}/wftask/"
            f"?{task_id=}"
        ),
        json=workflow_task,
    )
    workflow_task = check_response(res, expected_status_code=201)

    if batch:
        return PrintInterface(retcode=0, data=str(workflow_task["id"]))
    else:
        return RichJsonInterface(retcode=0, data=workflow_task)


def patch_workflowtask(
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

    res = client.patch(
        (
            f"{settings.BASE_URL}/project/{project_id}/"
            f"workflow/{workflow_id}/wftask/{workflow_task_id}/"
        ),
        json=payload,
    )
    workflow_task = check_response(res, expected_status_code=200)

    return RichJsonInterface(retcode=0, data=workflow_task)


def delete_workflowtask(
    client: AuthClient,
    *,
    project_id: int,
    workflow_id: int,
    workflow_task_id: int,
) -> BaseInterface:

    res = client.delete(
        f"{settings.BASE_URL}/project/{project_id}/"
        f"workflow/{workflow_id}/wftask/{workflow_task_id}/"
    )
    check_response(res, expected_status_code=204)
    return PrintInterface(retcode=0, data="")


def patch_workflow(
    client: AuthClient,
    *,
    project_id: int,
    workflow_id: int,
    new_name: str,
) -> BaseInterface:

    workflow_update = dict(name=new_name)

    res = client.patch(
        f"{settings.BASE_URL}/project/{project_id}/workflow/{workflow_id}/",
        json=workflow_update,
    )
    new_workflow = check_response(res, expected_status_code=200)
    return RichJsonInterface(retcode=0, data=new_workflow)


def workflow_apply(
    client: AuthClient,
    *,
    project_id: int,
    workflow_id: int,
    input_dataset_id: int,
    output_dataset_id: int,
    first_task_index: Optional[int] = None,
    last_task_index: Optional[int] = None,
    worker_init: Optional[str] = None,
    batch: bool = False,
) -> BaseInterface:
    apply_wf_create = dict(
        workflow_id=workflow_id,
        input_dataset_id=input_dataset_id,
        output_dataset_id=output_dataset_id,
    )
    # Prepare ApplyWorkflowCreate object, without None attributes
    if worker_init is not None:
        apply_wf_create["worker_init"] = worker_init
    if first_task_index is not None:
        apply_wf_create["first_task_index"] = first_task_index
    if last_task_index is not None:
        apply_wf_create["last_task_index"] = last_task_index

    # Prepare query parameters
    query_parameters = (
        f"input_dataset_id={input_dataset_id}"
        f"&output_dataset_id={output_dataset_id}"
    )

    res = client.post(
        (
            f"{settings.BASE_URL}/project/{project_id}/workflow/{workflow_id}/"
            f"apply/?{query_parameters}"
        ),
        json=apply_wf_create,
    )
    apply_wf_read = check_response(res, expected_status_code=202)

    if batch:
        return PrintInterface(retcode=0, data=apply_wf_read["id"])
    else:
        return RichJsonInterface(retcode=0, data=apply_wf_read)


def workflow_import(
    client: AuthClient, *, project_id: int, json_file: str, batch: bool = False
) -> BaseInterface:
    with Path(json_file).open("r") as f:
        workflow = json.load(f)

    res = client.post(
        f"{settings.BASE_URL}/project/{project_id}/workflow/import/",
        json=workflow,
    )
    wf_read = check_response(res, expected_status_code=201)

    warnings = [
        workflow_task["task"]["source"]
        for workflow_task in wf_read["task_list"]
        if workflow_task["task"]["owner"]
    ]
    if warnings:
        sources_str = ", ".join([f"'{s}'" for s in warnings])
        logging.warning(
            "This workflow includes custom tasks (the ones with sources: "
            f"{sources_str}), which are not meant to be portable; "
            "importing this workflow may not work as expected."
        )

    if batch:
        datastr = f"{wf_read['id']}"
        for wftask in wf_read["task_list"]:
            datastr += f" {wftask['id']}"
        return PrintInterface(retcode=0, data=datastr)
    else:
        return RichJsonInterface(retcode=0, data=wf_read)


def workflow_export(
    client: AuthClient,
    *,
    project_id: int,
    workflow_id: int,
    json_file: str,
) -> BaseInterface:
    res = client.get(
        (
            f"{settings.BASE_URL}/project/{project_id}/"
            f"workflow/{workflow_id}/export/"
        ),
    )
    workflow = check_response(res, expected_status_code=200)

    warnings = [
        workflow_task["task"]["source"]
        for workflow_task in workflow["task_list"]
        if not workflow_task["task"]["source"].startswith(
            ("pip_local:", "pip_remote:")
        )
    ]
    if warnings:
        sources_str = ", ".join([f"'{s}'" for s in warnings])
        logging.warning(
            "This workflow includes custom tasks (the ones with sources: "
            f"{sources_str}), which are not meant to be portable; "
            "re-importing this workflow may not work as expected."
        )

    with Path(json_file).open("w") as f:
        json.dump(workflow, f, indent=2)
    return PrintInterface(
        retcode=0, data=f"Workflow {workflow_id} exported at {json_file}"
    )
