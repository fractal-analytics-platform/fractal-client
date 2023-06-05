from typing import Optional

from httpx import AsyncClient

from ..authclient import AuthClient
from ..config import __VERSION__
from ..config import settings
from ..interface import BaseInterface
from ..interface import PrintInterface
from ._dataset import delete_dataset
from ._dataset import delete_resource
from ._dataset import get_dataset
from ._dataset import patch_dataset
from ._dataset import post_dataset
from ._dataset import post_resource
from ._job import get_job
from ._job import get_job_list
from ._job import get_job_logs
from ._job import stop_job
from ._project import delete_project
from ._project import get_project
from ._project import get_project_list
from ._project import patch_project
from ._project import post_project
from ._task import delete_task
from ._task import get_task_list
from ._task import patch_task
from ._task import post_task
from ._task import task_collect_pip
from ._task import task_collection_check
from ._user import user_delete
from ._user import user_edit
from ._user import user_list
from ._user import user_register
from ._user import user_show
from ._user import user_whoami
from ._workflow import delete_workflow
from ._workflow import delete_workflowtask
from ._workflow import get_workflow
from ._workflow import get_workflow_list
from ._workflow import patch_workflow
from ._workflow import patch_workflowtask
from ._workflow import post_workflow
from ._workflow import post_workflowtask
from ._workflow import workflow_apply
from ._workflow import workflow_export
from ._workflow import workflow_import


class NoCommandError(ValueError):
    pass


async def project(
    client: AuthClient,
    subcmd: str,
    name: Optional[str] = None,
    dataset: Optional[str] = None,
    project_id: Optional[int] = None,
    dataset_name: Optional[str] = None,
    metadata: Optional[str] = None,
    type: Optional[str] = None,
    new_name: Optional[str] = None,
    make_read_only: bool = False,
    remove_read_only: bool = False,
    verbose: bool = False,
    batch: bool = False,
) -> BaseInterface:
    if subcmd == "new":
        iface = await post_project(
            client, name=name, dataset=dataset, batch=batch
        )
    elif subcmd == "show":
        iface = await get_project(client, project_id=project_id)
    elif subcmd == "list":
        iface = await get_project_list(client)
    elif subcmd == "edit":
        iface = await patch_project(
            client,
            project_id=project_id,
            new_name=new_name,
            make_read_only=make_read_only,
            remove_read_only=remove_read_only,
        )
    elif subcmd == "add-dataset":
        iface = await post_dataset(
            client,
            project_id=project_id,
            dataset_name=dataset_name,
            metadata_filename=metadata,
            type=type,
            batch=batch,
        )
    elif subcmd == "delete":
        iface = await delete_project(client, project_id=project_id)
    else:
        raise NoCommandError(f"Command project {subcmd} not found")

    return iface


async def dataset(
    client: AuthClient,
    subcmd: str,
    project_id: int,
    dataset_id: int,
    resource_id: Optional[int] = None,
    path: Optional[str] = None,
    new_name: Optional[str] = None,
    new_type: Optional[str] = None,
    meta_file: Optional[str] = None,
    make_read_only: bool = False,
    remove_read_only: bool = False,
    verbose: bool = False,
    batch: bool = False,
) -> BaseInterface:
    if subcmd == "show":
        iface = await get_dataset(
            client, project_id=project_id, dataset_id=dataset_id
        )
    elif subcmd == "add-resource":
        iface = await post_resource(
            client,
            project_id=project_id,
            dataset_id=dataset_id,
            path=path,
            batch=batch,
        )
    elif subcmd == "rm-resource":
        iface = await delete_resource(
            client,
            project_id=project_id,
            dataset_id=dataset_id,
            resource_id=resource_id,
        )
    elif subcmd == "edit":
        iface = await patch_dataset(
            client,
            project_id=project_id,
            dataset_id=dataset_id,
            new_name=new_name,
            new_type=new_type,
            meta_file=meta_file,
            make_read_only=make_read_only,
            remove_read_only=remove_read_only,
        )
    elif subcmd == "delete":
        iface = await delete_dataset(
            client, project_id=project_id, dataset_id=dataset_id
        )
    else:
        raise NoCommandError(f"Command dataset {subcmd} not found")
    return iface


async def task(
    client: AuthClient,
    subcmd: str,
    package: Optional[str] = None,
    python_version: Optional[str] = None,
    package_version: Optional[str] = None,
    package_extras: Optional[str] = None,
    state_id: Optional[str] = None,
    name: Optional[str] = None,
    command: Optional[str] = None,
    source: Optional[str] = None,
    input_type: Optional[str] = None,
    output_type: Optional[str] = None,
    version: Optional[str] = None,
    default_args_file: Optional[str] = None,
    meta_file: Optional[str] = None,
    task_id_or_name: Optional[str] = None,
    new_name: Optional[str] = None,
    new_command: Optional[str] = None,
    new_input_type: Optional[str] = None,
    new_output_type: Optional[str] = None,
    new_version: Optional[str] = None,
    include_logs: bool = False,
    do_not_separate_logs: bool = False,
    batch: bool = False,
    verbose: bool = False,
) -> BaseInterface:

    if subcmd == "list":
        iface = await get_task_list(client)
    elif subcmd == "collect":
        iface = await task_collect_pip(
            client,
            package=package,
            package_version=package_version,
            python_version=python_version,
            package_extras=package_extras,
            batch=batch,
        )
    elif subcmd == "check-collection":
        iface = await task_collection_check(
            client,
            state_id=state_id,
            include_logs=include_logs,
            do_not_separate_logs=do_not_separate_logs,
        )
    elif subcmd == "new":
        iface = await post_task(
            client,
            name=name,
            command=command,
            source=source,
            input_type=input_type,
            output_type=output_type,
            version=version,
            default_args_file=default_args_file,
            meta_file=meta_file,
            batch=batch,
        )
    elif subcmd == "edit":
        iface = await patch_task(
            client,
            task_id_or_name=task_id_or_name,
            version=version,
            new_name=new_name,
            new_command=new_command,
            new_input_type=new_input_type,
            new_output_type=new_output_type,
            new_version=new_version,
            default_args_file=default_args_file,
            meta_file=meta_file,
        )
    elif subcmd == "delete":
        iface = await delete_task(client, task_id=task_id_or_name)
    else:
        raise NoCommandError(f"Command task {subcmd} not found")
    return iface


async def workflow(
    client: AuthClient,
    subcmd: str,
    batch: bool = False,
    verbose: bool = False,
    **kwargs,
) -> BaseInterface:
    if subcmd == "show":
        iface = await get_workflow(client, **kwargs)
    elif subcmd == "new":
        iface = await post_workflow(client, batch=batch, **kwargs)
    elif subcmd == "list":
        iface = await get_workflow_list(client, batch=batch, **kwargs)
    elif subcmd == "edit":
        iface = await patch_workflow(client, **kwargs)
    elif subcmd == "delete":
        iface = await delete_workflow(client, **kwargs)
    elif subcmd == "add-task":
        iface = await post_workflowtask(client, batch=batch, **kwargs)
    elif subcmd == "edit-task":
        iface = await patch_workflowtask(client, **kwargs)
    elif subcmd == "rm-task":
        iface = await delete_workflowtask(client, **kwargs)
    elif subcmd == "apply":
        iface = await workflow_apply(client, **kwargs)
    elif subcmd == "import":
        iface = await workflow_import(client, batch=batch, **kwargs)
    elif subcmd == "export":
        iface = await workflow_export(client, **kwargs)
    else:
        raise NoCommandError(f"Command workflow {subcmd} not found")
    return iface


async def job(
    client: AuthClient, subcmd: str, batch: bool = False, **kwargs
) -> BaseInterface:
    if subcmd == "list":
        iface = await get_job_list(client, batch=batch, **kwargs)
    elif subcmd == "show":
        iface = await get_job(client, batch=batch, **kwargs)
    elif subcmd == "download-logs":
        iface = await get_job_logs(client, **kwargs)
    elif subcmd == "stop":
        iface = await stop_job(client, batch=batch, **kwargs)
    else:
        raise NoCommandError(f"Command job {subcmd} not found")
    return iface


async def version(client: AsyncClient, **kwargs) -> PrintInterface:
    res = await client.get(f"{settings.FRACTAL_SERVER}/api/alive/")
    data = res.json()

    return PrintInterface(
        retcode=0,
        data=(
            f"Fractal client\n\tversion: {__VERSION__}\n"
            "Fractal server:\n"
            f"\turl: {settings.FRACTAL_SERVER}"
            f"\tdeployment type: {data['deployment_type']}"
            f"\tversion: {data['version']}"
        ),
    )


async def user(
    client: AuthClient, subcmd: str, batch: bool = False, **kwargs
) -> BaseInterface:
    if subcmd == "register":
        iface = await user_register(client, batch=batch, **kwargs)
    elif subcmd == "list":
        iface = await user_list(client, **kwargs)
    elif subcmd == "show":
        iface = await user_show(client, **kwargs)
    elif subcmd == "edit":
        iface = await user_edit(client, **kwargs)
    elif subcmd == "delete":
        iface = await user_delete(client, **kwargs)
    elif subcmd == "whoami":
        iface = await user_whoami(client, batch=batch, **kwargs)
    else:
        raise NoCommandError(f"Command user {subcmd} not found")

    return iface
