from httpx import AsyncClient

from ..authclient import AuthClient
from ..config import __VERSION__
from ..config import settings
from ..interface import BaseInterface
from ..interface import PrintInterface
from ..interface import RichJsonInterface
from ..response import check_response
from ._dataset import dataset_add_resource
from ._dataset import dataset_delete_resource
from ._dataset import dataset_edit
from ._dataset import dataset_show
from ._job import job_download_logs
from ._job import job_list
from ._job import job_status
from ._project import project_add_dataset
from ._project import project_create
from ._project import project_list
from ._project import project_show
from ._task import task_collect_pip
from ._task import task_collection_check
from ._task import task_edit
from ._task import task_list
from ._workflow import workflow_add_task
from ._workflow import workflow_apply
from ._workflow import workflow_delete
from ._workflow import workflow_edit
from ._workflow import workflow_edit_task
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
    elif subcmd == "add-dataset":
        iface = await project_add_dataset(
            client,
            batch=batch,
            project_id=kwargs.pop("project_id"),
            dataset_name=kwargs.pop("dataset_name"),
            metadata_filename=kwargs.pop("metadata"),
            **kwargs,
        )

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
            glob_pattern=kwargs.pop("glob_pattern"),
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

        dataset_update_dict = kwargs
        iface = await dataset_edit(
            client,
            project_id=project_id,
            dataset_id=dataset_id,
            dataset_update_dict=dataset_update_dict,
        )
    return iface


async def register(
    client: AsyncClient,
    email: str,
    slurm_user: str,
    password: str = None,
    **kwargs,
) -> BaseInterface:
    from getpass import getpass

    if password is None:
        password = getpass()
        confirm_password = getpass("Confirm password: ")
        if password == confirm_password:
            res = await client.post(
                f"{settings.FRACTAL_SERVER}/auth/register",
                json=dict(
                    email=email, slurm_user=slurm_user, password=password
                ),
            )

            data = check_response(res, expected_status_code=201)
            iface = RichJsonInterface(retcode=0, data=data)
        else:
            iface = PrintInterface(retcode=1, data="Passwords do not match.")
        return iface

    res = await client.post(
        f"{settings.FRACTAL_SERVER}/auth/register",
        json=dict(email=email, slurm_user=slurm_user, password=password),
    )
    data = check_response(res, expected_status_code=201)
    iface = RichJsonInterface(retcode=0, data=data)

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
    elif subcmd == "edit":
        iface = await task_edit(client, **kwargs)
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
        iface = await workflow_add_task(client, **kwargs)
    elif subcmd == "edit-task":
        iface = await workflow_edit_task(client, **kwargs)
    elif subcmd == "rm-task":
        iface = await workflow_remove_task(client, **kwargs)
    elif subcmd == "apply":
        iface = await workflow_apply(client, **kwargs)
    return iface


async def job(
    client: AuthClient, subcmd: str, batch: bool = False, **kwargs
) -> BaseInterface:
    if subcmd == "list":
        iface = await job_list(client, batch=batch, **kwargs)
    elif subcmd == "status":
        iface = await job_status(client, batch=batch, **kwargs)
    elif subcmd == "download-logs":
        iface = await job_download_logs(client, **kwargs)
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
