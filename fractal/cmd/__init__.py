from httpx import Client

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


def project(
    client: AuthClient, subcmd: str, batch: bool = False, **kwargs
) -> BaseInterface:
    if subcmd == "new":
        iface = post_project(client, batch=batch, **kwargs)
    elif subcmd == "show":
        iface = get_project(client, **kwargs)
    elif subcmd == "list":
        iface = get_project_list(client, **kwargs)
    elif subcmd == "edit":
        iface = patch_project(client, **kwargs)
    elif subcmd == "add-dataset":
        iface = post_dataset(
            client,
            batch=batch,
            project_id=kwargs.pop("project_id"),
            dataset_name=kwargs.pop("dataset_name"),
            metadata_filename=kwargs.pop("metadata"),
            type=kwargs.pop("type"),
            **kwargs,
        )
    elif subcmd == "delete":
        iface = delete_project(client, **kwargs)
    else:
        raise NoCommandError(f"Command project {subcmd} not found")

    return iface


def dataset(
    client: AuthClient, subcmd: str, batch: bool = False, **kwargs
) -> BaseInterface:
    if subcmd == "show":
        iface = get_dataset(client, **kwargs)
    elif subcmd == "add-resource":
        iface = post_resource(
            client,
            batch=batch,
            project_id=kwargs.pop("project_id"),
            dataset_id=kwargs.pop("dataset_id"),
            path=kwargs.pop("path"),
            **kwargs,
        )
    elif subcmd == "rm-resource":
        iface = delete_resource(
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
        iface = patch_dataset(
            client,
            project_id=project_id,
            dataset_id=dataset_id,
            **kwargs,
        )
    elif subcmd == "delete":
        iface = delete_dataset(client, **kwargs)
    else:
        raise NoCommandError(f"Command dataset {subcmd} not found")
    return iface


def task(
    client: AuthClient, subcmd: str, batch: bool = False, **kwargs
) -> BaseInterface:
    if subcmd == "list":
        iface = get_task_list(client)
    elif subcmd == "collect":
        iface = task_collect_pip(client, batch=batch, **kwargs)
    elif subcmd == "check-collection":
        iface = task_collection_check(client, **kwargs)
    elif subcmd == "new":
        iface = post_task(client, batch=batch, **kwargs)
    elif subcmd == "edit":
        iface = patch_task(client, **kwargs)
    elif subcmd == "delete":
        iface = delete_task(client, **kwargs)
    else:
        raise NoCommandError(f"Command task {subcmd} not found")
    return iface


def workflow(
    client: AuthClient, subcmd: str, batch: bool = False, **kwargs
) -> BaseInterface:
    if subcmd == "show":
        iface = get_workflow(client, **kwargs)
    elif subcmd == "new":
        iface = post_workflow(client, batch=batch, **kwargs)
    elif subcmd == "list":
        iface = get_workflow_list(client, batch=batch, **kwargs)
    elif subcmd == "edit":
        iface = patch_workflow(client, **kwargs)
    elif subcmd == "delete":
        iface = delete_workflow(client, **kwargs)
    elif subcmd == "add-task":
        iface = post_workflowtask(client, batch=batch, **kwargs)
    elif subcmd == "edit-task":
        iface = patch_workflowtask(client, **kwargs)
    elif subcmd == "rm-task":
        iface = delete_workflowtask(client, **kwargs)
    elif subcmd == "apply":
        iface = workflow_apply(client, **kwargs)
    elif subcmd == "import":
        iface = workflow_import(client, batch=batch, **kwargs)
    elif subcmd == "export":
        iface = workflow_export(client, **kwargs)
    else:
        raise NoCommandError(f"Command workflow {subcmd} not found")
    return iface


def job(
    client: AuthClient, subcmd: str, batch: bool = False, **kwargs
) -> BaseInterface:
    if subcmd == "list":
        iface = get_job_list(client, batch=batch, **kwargs)
    elif subcmd == "show":
        iface = get_job(client, batch=batch, **kwargs)
    elif subcmd == "download-logs":
        iface = get_job_logs(client, **kwargs)
    elif subcmd == "stop":
        iface = stop_job(client, batch=batch, **kwargs)
    else:
        raise NoCommandError(f"Command job {subcmd} not found")
    return iface


def version(client: Client, **kwargs) -> PrintInterface:
    res = client.get(f"{settings.FRACTAL_SERVER}/api/alive/")
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


def user(
    client: AuthClient, subcmd: str, batch: bool = False, **kwargs
) -> BaseInterface:
    if subcmd == "register":
        iface = user_register(client, batch=batch, **kwargs)
    elif subcmd == "list":
        iface = user_list(client, **kwargs)
    elif subcmd == "show":
        iface = user_show(client, **kwargs)
    elif subcmd == "edit":
        iface = user_edit(client, **kwargs)
    elif subcmd == "delete":
        iface = user_delete(client, **kwargs)
    elif subcmd == "whoami":
        iface = user_whoami(client, batch=batch, **kwargs)
    else:
        raise NoCommandError(f"Command user {subcmd} not found")

    return iface
