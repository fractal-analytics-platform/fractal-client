from ..authclient import AuthClient
from ..config import __VERSION__
from ..config import settings
from ..interface import BaseInterface
from ..interface import PrintInterface
from ._dataset import dataset_add_resource
from ._dataset import dataset_edit
from ._dataset import dataset_show
from ._project import project_add_dataset
from ._project import project_create
from ._project import project_list
from ._project import project_show
from ._task import task_list
from ._task import task_new


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
            project_id=kwargs.pop("project_id"),
            dataset_id=kwargs.pop("dataset_id"),
            path=kwargs.pop("path"),
            glob_pattern=kwargs.pop("glob_pattern"),
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


async def register():
    raise NotImplementedError


async def task(
    client: AuthClient, subcmd: str, batch: bool = False, **kwargs
) -> BaseInterface:
    iface = None
    if subcmd == "list":
        iface = await task_list(client, **kwargs)
    elif subcmd == "new":
        iface = await task_new(client, **kwargs)
        pass
    elif subcmd == "edit":
        pass
    elif subcmd == "add-subtask":
        pass
    elif subcmd == "apply":
        pass
    return iface


async def version(client: AuthClient, **kwargs) -> PrintInterface:
    res = await client.get(f"{settings.FRACTAL_SERVER}/api/alive/")
    data = res.json()

    return PrintInterface(
        retcode=0,
        data=(
            f"Fractal client\n\tversion: {__VERSION__}\n"
            "Fractal server:"
            f"\turl: {settings.FRACTAL_SERVER}"
            f"\tdeployment type: {data['deployment_type']}"
            f"\tversion: {data['version']}"
        ),
    )
