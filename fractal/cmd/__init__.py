from httpx import AsyncClient

from ..authclient import AuthClient
from ..config import __VERSION__
from ..config import settings
from ..interface import BaseInterface
from ..interface import PrintInterface
from ._dataset import dataset_add_resource
from ._dataset import dataset_delete
from ._dataset import dataset_delete_resource
from ._dataset import dataset_edit
from ._dataset import dataset_show
from ._job import job_download_logs
from ._job import job_list
from ._job import job_show
from ._project import project_add_dataset
from ._project import project_create
from ._project import project_delete
from ._project import project_edit
from ._project import project_list
from ._project import project_show
from ._task import task_collect_pip
from ._task import task_collection_check
from ._task import task_delete
from ._task import task_edit
from ._task import task_list
from ._task import task_new
from ._user import user_delete
from ._user import user_edit
from ._user import user_list
from ._user import user_register
from ._user import user_show
from ._user import user_whoami
from ._workflow import workflow_add_task
from ._workflow import workflow_apply
from ._workflow import workflow_delete
from ._workflow import workflow_edit
from ._workflow import workflow_edit_task
from ._workflow import workflow_export
from ._workflow import workflow_import
from ._workflow import workflow_list
from ._workflow import workflow_new
from ._workflow import workflow_remove_task
from ._workflow import workflow_show


class NoCommandError(ValueError):
    pass


async def project(
    client: AuthClient, subcmd: str, batch: bool = False, **kwargs
) -> BaseInterface:
    if subcmd == "new":
        iface = await project_create(client, batch=batch, **kwargs)
    elif subcmd == "show":
        iface = await project_show(client, **kwargs)
    elif subcmd == "list":
        iface = await project_list(client, **kwargs)
    elif subcmd == "edit":
        iface = await project_edit(client, **kwargs)
    elif subcmd == "add-dataset":
        iface = await project_add_dataset(
            client,
            batch=batch,
            project_id=kwargs.pop("project_id"),
            dataset_name=kwargs.pop("dataset_name"),
            metadata_filename=kwargs.pop("metadata"),
            **kwargs,
        )
    elif subcmd == "delete":
        iface = await project_delete(client, **kwargs)
    else:
        raise NoCommandError(f"Command project {subcmd} not found")

    return iface


async def dataset(
    client: AuthClient, subcmd: str, batch: bool = False, **kwargs
) -> BaseInterface:
    if subcmd == "show":
        iface = await dataset_show(client, **kwargs)
    elif subcmd == "add-resource":
        iface = await dataset_add_resource(
            client,
            batch=batch,
            project_id=kwargs.pop("project_id"),
            dataset_id=kwargs.pop("dataset_id"),
            path=kwargs.pop("path"),
            **kwargs,
        )
    elif subcmd == "rm-resource":
        iface = await dataset_delete_resource(
            client,
            batch=batch,
            project_id=kwargs.pop("project_id"),
            dataset_id=kwargs.pop("dataset_id"),
            resource_id=kwargs.pop("resource_id"),
            **kwargs,
        )
    elif subcmd == "edit":
        project_id = int(kwargs.pop("project_id"))
        dataset_id = int(kwargs.pop("dataset_id"))
        iface = await dataset_edit(
            client,
            project_id=project_id,
            dataset_id=dataset_id,
            **kwargs,
        )
    elif subcmd == "delete":
        iface = await dataset_delete(client, **kwargs)
    else:
        raise NoCommandError(f"Command dataset {subcmd} not found")
    return iface


async def task(
    client: AuthClient, subcmd: str, batch: bool = False, **kwargs
) -> BaseInterface:
    if subcmd == "list":
        iface = await task_list(client, **kwargs)
    elif subcmd == "collect":
        iface = await task_collect_pip(client, batch=batch, **kwargs)
    elif subcmd == "check-collection":
        iface = await task_collection_check(client, **kwargs)
    elif subcmd == "new":
        iface = await task_new(client, batch=batch, **kwargs)
    elif subcmd == "edit":
        iface = await task_edit(client, **kwargs)
    elif subcmd == "delete":
        iface = await task_delete(client, **kwargs)
    else:
        raise NoCommandError(f"Command task {subcmd} not found")
    return iface


async def workflow(
    client: AuthClient, subcmd: str, batch: bool = False, **kwargs
) -> BaseInterface:
    if subcmd == "show":
        iface = await workflow_show(client, **kwargs)
    elif subcmd == "new":
        iface = await workflow_new(client, batch=batch, **kwargs)
    elif subcmd == "list":
        iface = await workflow_list(client, batch=batch, **kwargs)
    elif subcmd == "edit":
        iface = await workflow_edit(client, **kwargs)
    elif subcmd == "delete":
        iface = await workflow_delete(client, **kwargs)
    elif subcmd == "add-task":
        iface = await workflow_add_task(client, batch=batch, **kwargs)
    elif subcmd == "edit-task":
        iface = await workflow_edit_task(client, **kwargs)
    elif subcmd == "rm-task":
        iface = await workflow_remove_task(client, **kwargs)
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
        iface = await job_list(client, batch=batch, **kwargs)
    elif subcmd == "show":
        iface = await job_show(client, batch=batch, **kwargs)
    elif subcmd == "download-logs":
        iface = await job_download_logs(client, **kwargs)
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
