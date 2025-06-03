import json
import logging
import sys
from pathlib import Path

from ..authclient import AuthClient
from ..interface import Interface
from ..response import check_response
from ._aux_task_caching import FractalCacheError
from ._aux_task_caching import get_task_id_from_cache


def post_workflow(
    client: AuthClient, *, name: str, project_id: int, batch: bool = False
) -> Interface:
    workflow = dict(
        name=name,
    )
    res = client.post(
        f"api/v2/project/{project_id}/workflow/",
        json=workflow,
    )
    workflow = check_response(res, expected_status_code=201)
    if batch:
        return Interface(retcode=0, data=workflow["id"])
    else:
        return Interface(retcode=0, data=workflow)


def get_workflow_list(
    client: AuthClient, *, project_id: int, batch: bool = False
) -> Interface:

    res = client.get(f"api/v2/project/{project_id}/workflow/")
    workflow_list = check_response(res, expected_status_code=200)
    return Interface(retcode=0, data=workflow_list)


def delete_workflow(
    client: AuthClient, *, project_id: int, workflow_id: int
) -> Interface:
    res = client.delete(f"api/v2/project/{project_id}/workflow/{workflow_id}/")
    check_response(res, expected_status_code=204)
    return Interface(retcode=0, data="")


def get_workflow(
    client: AuthClient, *, project_id: int, workflow_id: int
) -> Interface:
    res = client.get(f"api/v2/project/{project_id}/workflow/{workflow_id}/")
    workflow = check_response(res, expected_status_code=200)
    return Interface(retcode=0, data=workflow)


def post_workflowtask(
    client: AuthClient,
    *,
    project_id: int,
    workflow_id: int,
    type_filters: str | None = None,
    args_non_parallel: str | None = None,
    args_parallel: str | None = None,
    meta_non_parallel: str | None = None,
    meta_parallel: str | None = None,
    task_id: int | None = None,
    task_name: str | None = None,
    task_version: str | None = None,
    batch: bool = False,
    order: int | None = None,
) -> Interface:

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

    if type_filters is not None:
        with Path(type_filters).open("r") as f:
            type_filters = json.load(f)
            workflow_task["type_filters"] = type_filters

    if args_non_parallel is not None:
        with Path(args_non_parallel).open("r") as f:
            a_n_p = json.load(f)
            workflow_task["args_non_parallel"] = a_n_p

    if args_parallel is not None:
        with Path(args_parallel).open("r") as f:
            a_p = json.load(f)
            workflow_task["args_parallel"] = a_p

    if meta_non_parallel is not None:
        with Path(meta_non_parallel).open("r") as f:
            m_n_p = json.load(f)
            workflow_task["meta_non_parallel"] = m_n_p

    if meta_parallel is not None:
        with Path(meta_parallel).open("r") as f:
            m_p = json.load(f)
            workflow_task["meta_parallel"] = m_p

    res = client.post(
        (
            f"api/v2/project/{project_id}/"
            f"workflow/{workflow_id}/wftask/"
            f"?{task_id=}"
        ),
        json=workflow_task,
    )
    workflow_task = check_response(res, expected_status_code=201)

    if batch:
        return Interface(retcode=0, data=str(workflow_task["id"]))
    else:
        return Interface(retcode=0, data=workflow_task)


def patch_workflowtask(
    client: AuthClient,
    *,
    project_id: int,
    workflow_id: int,
    workflow_task_id: int,
    type_filters: str | None = None,
    args_non_parallel: str | None = None,
    args_parallel: str | None = None,
    meta_non_parallel: str | None = None,
    meta_parallel: str | None = None,
) -> Interface:

    payload = {}
    if type_filters is not None:
        with Path(type_filters).open("r") as f:
            filters = json.load(f)
            payload["type_filters"] = filters

    if args_non_parallel is not None:
        with Path(args_non_parallel).open("r") as f:
            a_n_p = json.load(f)
            payload["args_non_parallel"] = a_n_p

    if args_parallel is not None:
        with Path(args_parallel).open("r") as f:
            a_p = json.load(f)
            payload["args_parallel"] = a_p

    if meta_non_parallel is not None:
        with Path(meta_non_parallel).open("r") as f:
            m_n_p = json.load(f)
            payload["meta_non_parallel"] = m_n_p

    if meta_parallel is not None:
        with Path(meta_parallel).open("r") as f:
            m_p = json.load(f)
            payload["meta_parallel"] = m_p

    res = client.patch(
        (
            f"api/v2/project/{project_id}/"
            f"workflow/{workflow_id}/wftask/{workflow_task_id}/"
        ),
        json=payload,
    )
    workflow_task = check_response(res, expected_status_code=200)

    return Interface(retcode=0, data=workflow_task)


def delete_workflowtask(
    client: AuthClient,
    *,
    project_id: int,
    workflow_id: int,
    workflow_task_id: int,
) -> Interface:

    res = client.delete(
        f"api/v2/project/{project_id}/"
        f"workflow/{workflow_id}/wftask/{workflow_task_id}/"
    )
    check_response(res, expected_status_code=204)
    return Interface(retcode=0, data="")


def patch_workflow(
    client: AuthClient,
    *,
    project_id: int,
    workflow_id: int,
    new_name: str,
) -> Interface:

    workflow_update = dict(name=new_name)

    res = client.patch(
        f"api/v2/project/{project_id}/workflow/{workflow_id}/",
        json=workflow_update,
    )
    new_workflow = check_response(res, expected_status_code=200)
    return Interface(retcode=0, data=new_workflow)


def workflow_import(
    client: AuthClient,
    *,
    project_id: int,
    json_file: str,
    workflow_name: str | None = None,
    batch: bool = False,
) -> Interface:
    with Path(json_file).open("r") as f:
        workflow = json.load(f)

    if workflow_name is not None:
        workflow["name"] = workflow_name

    res = client.post(
        f"api/v2/project/{project_id}/workflow/import/",
        json=workflow,
    )
    wf_read = check_response(res, expected_status_code=201)

    if batch:
        datastr = f"{wf_read['id']}"
        for wftask in wf_read["task_list"]:
            datastr += f" {wftask['id']}"
        return Interface(retcode=0, data=datastr)
    else:
        return Interface(retcode=0, data=wf_read)


def workflow_export(
    client: AuthClient,
    *,
    project_id: int,
    workflow_id: int,
    json_file: str,
) -> Interface:
    res = client.get(
        (f"api/v2/project/{project_id}/" f"workflow/{workflow_id}/export/"),
    )
    workflow = check_response(res, expected_status_code=200)

    with Path(json_file).open("w") as f:
        json.dump(workflow, f, indent=2)
    return Interface(
        retcode=0, data=f"Workflow {workflow_id} exported at {json_file}"
    )
