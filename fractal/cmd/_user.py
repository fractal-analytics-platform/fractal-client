from ..authclient import AuthClient
from ..config import settings
from ..interface import PrintInterface

def user_register():
    pass


def user_show():
    pass


async def user_whoami(
    client: AuthClient, **kwargs
) -> PrintInterface:
    res = await client.get(f"{settings.FRACTAL_SERVER}/auth/users/me")
    data = res.json()

    return PrintInterface(
        retcode=0,
        data=data,
    )