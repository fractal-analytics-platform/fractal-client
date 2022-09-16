from ..authclient import AuthClient
from ..config import __VERSION__
from ..config import settings
from ..interface import BaseInterface
from ..interface import PrintInterface
from ._project import project_add_dataset
from ._project import project_create
from ._project import project_list


async def project(
    client: AuthClient, subcmd: str, batch: bool = False, **kwargs
) -> BaseInterface:
    if subcmd == "new":
        return await project_create(client, batch=batch, **kwargs)
    elif subcmd == "list":
        return await project_list(client, **kwargs)
    elif subcmd == "add-dataset":
        return await project_add_dataset(
            client,
            kwargs.pop("project_id"),
            kwargs.pop("dataset_name"),
            metadata_filename=kwargs.pop("metadata"),
            **kwargs,
        )


async def register():
    raise NotImplementedError


async def dataset():
    raise NotImplementedError


async def task():
    raise NotImplementedError


async def version(client: AuthClient, **kwargs) -> PrintInterface:
    res = await client.get(f"{settings.FRACTAL_SERVER}/api/alive/")
    data = res.json()

    return PrintInterface(
        retcode=0,
        output=(
            f"Fractal client\n\tversion: {__VERSION__}\n"
            "Fractal server:"
            f"\turl: {settings.FRACTAL_SERVER}"
            f"\tdeployment type: {data['deployment_type']}"
            f"\tversion: {data['version']}"
        ),
    )
