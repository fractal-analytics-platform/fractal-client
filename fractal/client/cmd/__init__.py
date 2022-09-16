from ..authclient import AuthClient
from ..config import __VERSION__
from ..config import settings
from ._project import project_create


async def project(
    client: AuthClient, subcmd: str, batch: bool = False, **kwargs
):
    if subcmd == "new":
        data = await project_create(client, **kwargs)
    # project_list()

    if batch:
        print(data["id"])
    else:
        print(data)


async def register():
    raise NotImplementedError


async def dataset():
    raise NotImplementedError


async def task():
    raise NotImplementedError


async def version(client: AuthClient, **kwargs):
    res = await client.get(f"{settings.FRACTAL_SERVER}/api/alive/")
    data = res.json()

    print(f"Fractal client\n\tversion: {__VERSION__}")
    print()
    print("Fractal server:")
    print(f"\turl: {settings.FRACTAL_SERVER}")
    print(f"\tdeployment type: {data['deployment_type']}")
    print(f"\tversion: {data['version']}")
