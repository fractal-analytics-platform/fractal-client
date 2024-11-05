from httpx import Client
from httpx import ConnectError

from ..authclient import AuthClient
from ..config import settings
from ..interface import Interface
from ._dataset import delete_dataset
from ._dataset import get_dataset
from ._dataset import patch_dataset
from ._dataset import post_dataset
from ._group import group_get
from ._group import group_list
from ._group import group_new
from ._group import group_update
from ._job import get_job
from ._job import get_job_list
from ._job import get_job_logs
from ._job import job_submit
from ._job import stop_job
from ._project import delete_project
from ._project import get_project
from ._project import get_project_list
from ._project import patch_project
from ._project import post_project
from ._task import get_task_list
from ._task import patch_task
from ._task import post_task
from ._task import task_collect_custom
from ._task import task_collect_pip
from ._task import task_collection_check
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
from ._workflow import workflow_export
from ._workflow import workflow_import
from fractal_client import __VERSION__


class NoCommandError(ValueError):
    pass


def get_kwargs(_parameters, _kwargs):
    return {k: _kwargs.get(k) for k in _parameters if k in _kwargs}


def project(
    client: AuthClient,
    subcmd: str,
    batch: bool = False,
    **kwargs,
) -> Interface:

    if subcmd == "new":
        parameters = ["name"]
        function_kwargs = get_kwargs(parameters, kwargs)
        iface = post_project(client, batch=batch, **function_kwargs)
    elif subcmd == "show":
        parameters = ["project_id"]
        function_kwargs = get_kwargs(parameters, kwargs)
        iface = get_project(client, **function_kwargs)
    elif subcmd == "list":
        iface = get_project_list(client)
    elif subcmd == "edit":
        parameters = [
            "project_id",
            "new_name",
        ]
        function_kwargs = get_kwargs(parameters, kwargs)
        iface = patch_project(client, **function_kwargs)
    elif subcmd == "add-dataset":
        parameters = ["project_id", "dataset_name", "zarr_dir", "filters"]
        function_kwargs = get_kwargs(parameters, kwargs)
        iface = post_dataset(client, batch=batch, **function_kwargs)
    elif subcmd == "delete":
        parameters = ["project_id"]
        function_kwargs = get_kwargs(parameters, kwargs)
        iface = delete_project(client, **function_kwargs)
    else:
        raise NoCommandError(f"Command 'project {subcmd}' not found")

    return iface


def dataset(
    client: AuthClient,
    subcmd: str,
    batch: bool = False,
    **kwargs,
) -> Interface:
    if subcmd == "show":
        parameters = ["project_id", "dataset_id"]
        function_kwargs = get_kwargs(parameters, kwargs)
        iface = get_dataset(client, **function_kwargs)
    elif subcmd == "edit":
        parameters = [
            "project_id",
            "dataset_id",
            "new_name",
            "filters",
        ]
        function_kwargs = get_kwargs(parameters, kwargs)
        iface = patch_dataset(client, **function_kwargs)
    elif subcmd == "delete":
        parameters = ["project_id", "dataset_id"]
        function_kwargs = get_kwargs(parameters, kwargs)
        iface = delete_dataset(client, **function_kwargs)
    else:
        raise NoCommandError(f"Command 'dataset {subcmd}' not found")
    return iface


def task(
    client: AuthClient,
    subcmd: str,
    batch: bool = False,
    **kwargs,
) -> Interface:

    if subcmd == "list":
        iface = get_task_list(client)
    elif subcmd == "collect":
        parameters = [
            "package",
            "package_version",
            "python_version",
            "package_extras",
            "pinned_dependency",
            "private",
        ]
        function_kwargs = get_kwargs(parameters, kwargs)
        iface = task_collect_pip(client, batch=batch, **function_kwargs)
    elif subcmd == "collect-custom":
        parameters = [
            "label",
            "python_interpreter",
            "manifest",
            "version",
            "package_name",
            "package_root",
            "private",
        ]
        function_kwargs = get_kwargs(parameters, kwargs)
        iface = task_collect_custom(client, batch=batch, **function_kwargs)
    elif subcmd == "check-collection":
        parameters = ["state_id", "include_logs", "do_not_separate_logs"]
        function_kwargs = get_kwargs(parameters, kwargs)
        iface = task_collection_check(client, **function_kwargs)
    elif subcmd == "new":
        parameters = [
            "name",
            "version",
            "command_non_parallel",
            "command_parallel",
            "meta_non_parallel",
            "meta_parallel",
            "args_schema_non_parallel",
            "args_schema_parallel",
            "args_schema_version",
            "private",
        ]
        function_kwargs = get_kwargs(parameters, kwargs)
        iface = post_task(client, batch=batch, **function_kwargs)
    elif subcmd == "edit":
        parameters = [
            "id",
            "name",
            "version",
            "command_non_parallel",
            "command_parallel",
            "input_types",
            "output_types",
        ]
        function_kwargs = get_kwargs(parameters, kwargs)
        iface = patch_task(client, **function_kwargs)
    else:
        raise NoCommandError(f"Command 'task {subcmd}' not found")
    return iface


def workflow(
    client: AuthClient,
    subcmd: str,
    batch: bool = False,
    **kwargs,
) -> Interface:
    if subcmd == "show":
        parameters = ["project_id", "workflow_id"]
        function_kwargs = get_kwargs(parameters, kwargs)
        iface = get_workflow(client, **function_kwargs)
    elif subcmd == "new":
        parameters = ["name", "project_id"]
        function_kwargs = get_kwargs(parameters, kwargs)
        iface = post_workflow(client, batch=batch, **function_kwargs)
    elif subcmd == "list":
        parameters = ["project_id"]
        function_kwargs = get_kwargs(parameters, kwargs)
        iface = get_workflow_list(client, batch=batch, **function_kwargs)
    elif subcmd == "edit":
        parameters = ["project_id", "workflow_id", "new_name"]
        function_kwargs = get_kwargs(parameters, kwargs)
        iface = patch_workflow(client, **function_kwargs)
    elif subcmd == "delete":
        parameters = ["project_id", "workflow_id"]
        function_kwargs = get_kwargs(parameters, kwargs)
        iface = delete_workflow(client, **function_kwargs)
    elif subcmd == "add-task":
        parameters = [
            "project_id",
            "workflow_id",
            "task_id",
            "task_name",
            "task_version",
            "order",
            "args_non_parallel",
            "args_parallel",
            "input_filters",
            "meta_parallel",
            "meta_non_parallel",
        ]
        function_kwargs = get_kwargs(parameters, kwargs)
        iface = post_workflowtask(client, batch=batch, **function_kwargs)
    elif subcmd == "edit-task":
        parameters = [
            "project_id",
            "workflow_id",
            "workflow_task_id",
            "input_filters",
            "args_non_parallel",
            "args_parallel",
            "meta_parallel",
            "meta_non_parallel",
        ]
        function_kwargs = get_kwargs(parameters, kwargs)
        iface = patch_workflowtask(client, **function_kwargs)
    elif subcmd == "rm-task":
        parameters = ["project_id", "workflow_id", "workflow_task_id"]
        function_kwargs = get_kwargs(parameters, kwargs)
        iface = delete_workflowtask(client, **function_kwargs)
    elif subcmd == "import":
        parameters = ["project_id", "json_file", "workflow_name"]
        function_kwargs = get_kwargs(parameters, kwargs)
        iface = workflow_import(client, batch=batch, **function_kwargs)
    elif subcmd == "export":
        parameters = ["project_id", "workflow_id", "json_file"]
        function_kwargs = get_kwargs(parameters, kwargs)
        iface = workflow_export(client, **function_kwargs)
    else:
        raise NoCommandError(f"Command 'workflow {subcmd}' not found")
    return iface


def job(
    client: AuthClient,
    subcmd: str,
    batch: bool = False,
    **kwargs,
) -> Interface:
    if subcmd == "list":
        parameters = ["project_id"]
        function_kwargs = get_kwargs(parameters, kwargs)
        iface = get_job_list(client, batch=batch, **function_kwargs)
    elif subcmd == "show":
        parameters = ["project_id", "job_id", "do_not_separate_logs"]
        function_kwargs = get_kwargs(parameters, kwargs)
        iface = get_job(client, batch=batch, **function_kwargs)
    elif subcmd == "download-logs":
        parameters = ["project_id", "job_id", "output_folder"]
        function_kwargs = get_kwargs(parameters, kwargs)
        iface = get_job_logs(client, **function_kwargs)
    elif subcmd == "stop":
        parameters = ["project_id", "job_id"]
        function_kwargs = get_kwargs(parameters, kwargs)
        iface = stop_job(client, **function_kwargs)
    elif subcmd == "submit":
        parameters = [
            "project_id",
            "workflow_id",
            "dataset_id",
            "worker_init",
            "first_task_index",
            "last_task_index",
        ]
        function_kwargs = get_kwargs(parameters, kwargs)
        iface = job_submit(client, batch=batch, **function_kwargs)
    else:
        raise NoCommandError(f"Command 'job {subcmd}' not found")
    return iface


def version(client: Client, **kwargs) -> Interface:
    try:
        res = client.get(f"{settings.FRACTAL_SERVER}/api/alive/")
        data = res.json()
        server_str = (
            f"\turl: {settings.FRACTAL_SERVER}\tversion: {data['version']}"
        )
    except ConnectError:
        server_str = f"\tConnection to '{settings.FRACTAL_SERVER}' refused"

    return Interface(
        retcode=0,
        data=(
            f"Fractal client\n\tversion: {__VERSION__}\n"
            f"Fractal server:\n{server_str}"
        ),
    )


def user(
    client: AuthClient, subcmd: str, batch: bool = False, **kwargs
) -> Interface:
    if subcmd == "register":
        parameters = [
            "new_email",
            "new_password",
            "cache_dir",
            "slurm_user",
            "project_dir",
            "username",
            "ssh_settings_json",
            "superuser",
        ]
        function_kwargs = get_kwargs(parameters, kwargs)
        iface = user_register(client, batch=batch, **function_kwargs)
    elif subcmd == "list":
        iface = user_list(client)
    elif subcmd == "show":
        parameters = ["user_id"]
        function_kwargs = get_kwargs(parameters, kwargs)
        iface = user_show(client, **function_kwargs)
    elif subcmd == "edit":
        parameters = [
            "user_id",
            "new_email",
            "new_password",
            "new_username",
            "new_slurm_user",
            "new_cache_dir",
            "new_project_dir",
            "new_ssh_settings_json",
            "make_superuser",
            "remove_superuser",
            "make_verified",
            "remove_verified",
        ]
        function_kwargs = get_kwargs(parameters, kwargs)
        iface = user_edit(client, **function_kwargs)
    elif subcmd == "whoami":
        parameters = ["viewer_paths"]
        function_kwargs = get_kwargs(parameters, kwargs)
        iface = user_whoami(client, **function_kwargs, batch=batch)
    else:
        raise NoCommandError(f"Command 'user {subcmd}' not found")

    return iface


def group(
    client: AuthClient, subcmd: str, batch: bool = False, **kwargs
) -> Interface:
    if subcmd == "list":
        parameters = ["user_ids"]
        function_kwargs = get_kwargs(parameters, kwargs)
        iface = group_list(client, batch=batch, **function_kwargs)
    elif subcmd == "get":
        parameters = ["group_id"]
        function_kwargs = get_kwargs(parameters, kwargs)
        iface = group_get(client, **function_kwargs)
    elif subcmd == "new":
        parameters = ["name", "viewer_paths"]
        function_kwargs = get_kwargs(parameters, kwargs)
        iface = group_new(client, batch=batch, **function_kwargs)
    elif subcmd == "update":
        parameters = ["group_id", "new_user_ids", "new_viewer_paths"]
        function_kwargs = get_kwargs(parameters, kwargs)
        iface = group_update(client, **function_kwargs)
    else:
        raise NoCommandError(f"Command 'group {subcmd}' not found")

    return iface
