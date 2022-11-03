import logging

from ...common.models import WorkflowCreate
from ...common.models import WorkflowRead
from ...common.models import WorkflowTaskCreate
from ...common.models import WorkflowUpdate
from ..authclient import AuthClient
from ..config import settings
from ..interface import BaseInterface
from ..interface import PrintInterface
from ..interface import RichJsonInterface
from ..response import check_response


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
    task_id: int,
    order: int = None,
    **kwargs,
) -> RichJsonInterface:
    workflow_task = WorkflowTaskCreate(task_id=task_id, order=order)
    res = await client.post(
        f"{settings.BASE_URL}/workflow/{id}/add-task/",
        json=workflow_task.dict(),
    )

    workflow_task = check_response(
        res, expected_status_code=201, coerce=WorkflowRead
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
