import json
from typing import Optional

from ..authclient import AuthClient
from ..config import settings
from ..interface import RichJsonInterface
from ..response import check_response
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
    project_name: str,
    input_dataset_name: str,
    output_dataset_name: str,
    workflow_name: str,
    overwrite_input: bool,
    **kwargs,
) -> RichJsonInterface:

    res = await client.get(f"{settings.BASE_URL}/project/")
    projects = res.json()
    # Find project
    try:
        project = [p for p in projects if p["name"] == project_name][0]
        project_id = project["id"]
    except IndexError as e:
        raise IndexError(f"Project {project_name} not found", str(e))
    # Get dataset list
    dataset_list = [
        project["dataset_list"]
        for project in projects
        if project["id"] == project_id
    ][0]
    # Find I/O dataset
    try:
        input_dataset = [
            ds for ds in dataset_list if ds["name"] == input_dataset_name
        ][0]
        input_dataset_id = input_dataset["id"]
    except IndexError as e:
        raise IndexError(f"Dataset {input_dataset_name} not found", str(e))
    try:
        output_dataset = [
            ds for ds in dataset_list if ds["name"] == output_dataset_name
        ][0]
        output_dataset_id = output_dataset["id"]
    except IndexError as e:
        raise IndexError(f"Dataset {output_dataset_name} not found", str(e))
    # Extract ID for task
    res = await client.get(f"{settings.BASE_URL}/task/")
    data = res.json()
    task_list = [TaskRead(**item) for item in data]
    # Find task id
    try:
        workflow = [t for t in task_list if t.name == workflow_name][0]
        workflow_id = workflow.id
    except IndexError as e:
        raise IndexError(f"Task {workflow_name} not found", str(e))

    # Apply workflow

    from fractal.common.models import ApplyWorkflow

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
    wf_submitted = check_response(res, expected_status_code=200)
    return RichJsonInterface(retcode=0, data=wf_submitted)
