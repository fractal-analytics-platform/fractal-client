import json
import logging
import sys
from pathlib import Path
from typing import Optional

from ..authclient import AuthClient
from ..config import settings
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
        f"{settings.BASE_URL}/project/{project_id}/workflow/",
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

    res = client.get(f"{settings.BASE_URL}/project/{project_id}/workflow/")
    workflow_list = check_response(res, expected_status_code=200)
    return Interface(retcode=0, data=workflow_list)


def delete_workflow(
    client: AuthClient, *, project_id: int, workflow_id: int
) -> Interface:
    res = client.delete(
        f"{settings.BASE_URL}/project/{project_id}/workflow/{workflow_id}/"
    )
    check_response(res, expected_status_code=204)
    return Interface(retcode=0, data="")


def get_workflow(
    client: AuthClient, *, project_id: int, workflow_id: int
) -> Interface:
    res = client.get(
        f"{settings.BASE_URL}/project/{project_id}/workflow/{workflow_id}/"
    )
    workflow = check_response(res, expected_status_code=200)
    return Interface(retcode=0, data=workflow)


def post_workflowtask(
    client: AuthClient,
    *,
    project_id: int,
    workflow_id: int,
    input_filters: Optional[str] = None,
    args_non_parallel: Optional[str] = None,
    args_parallel: Optional[str] = None,
    meta_non_parallel: Optional[str] = None,
    meta_parallel: Optional[str] = None,
    task_id: Optional[int] = None,
    task_name: Optional[str] = None,
    task_version: Optional[str] = None,
    batch: bool = False,
    order: Optional[int] = None,
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

    if input_filters:
        with Path(input_filters).open("r") as f:
            i_filters = json.load(f)
            workflow_task["input_filters"] = i_filters

    if args_non_parallel:
        with Path(args_non_parallel).open("r") as f:
            a_n_p = json.load(f)
            workflow_task["args_non_parallel"] = a_n_p

    if args_parallel:
        with Path(args_parallel).open("r") as f:
            a_p = json.load(f)
            workflow_task["args_parallel"] = a_p

    if meta_non_parallel:
        with Path(meta_non_parallel).open("r") as f:
            m_n_p = json.load(f)
            workflow_task["meta_non_parallel"] = m_n_p

    if meta_parallel:
        with Path(meta_parallel).open("r") as f:
            m_p = json.load(f)
            workflow_task["meta_parallel"] = m_p

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
        return Interface(retcode=0, data=str(workflow_task["id"]))
    else:
        return Interface(retcode=0, data=workflow_task)


def patch_workflowtask(
    client: AuthClient,
    *,
    project_id: int,
    workflow_id: int,
    workflow_task_id: int,
    input_filters: Optional[str] = None,
    args_non_parallel: Optional[str] = None,
    args_parallel: Optional[str] = None,
    meta_non_parallel: Optional[str] = None,
    meta_parallel: Optional[str] = None,
) -> Interface:

    payload = {}
    if input_filters:
        with Path(input_filters).open("r") as f:
            input_filters = json.load(f)
            payload["input_filters"] = input_filters

    if args_non_parallel:
        with Path(args_non_parallel).open("r") as f:
            a_n_p = json.load(f)
            payload["args_non_parallel"] = a_n_p

    if args_parallel:
        with Path(args_parallel).open("r") as f:
            a_p = json.load(f)
            payload["args_parallel"] = a_p

    if meta_non_parallel:
        with Path(meta_non_parallel).open("r") as f:
            m_n_p = json.load(f)
            payload["meta_non_parallel"] = m_n_p

    if meta_parallel:
        with Path(meta_parallel).open("r") as f:
            m_p = json.load(f)
            payload["meta_parallel"] = m_p

    res = client.patch(
        (
            f"{settings.BASE_URL}/project/{project_id}/"
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
        f"{settings.BASE_URL}/project/{project_id}/"
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
        f"{settings.BASE_URL}/project/{project_id}/workflow/{workflow_id}/",
        json=workflow_update,
    )
    new_workflow = check_response(res, expected_status_code=200)
    return Interface(retcode=0, data=new_workflow)


def workflow_import() -> None:
    raise NotImplementedError


def workflow_export() -> None:
    raise NotImplementedError
