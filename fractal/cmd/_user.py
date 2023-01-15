from ..authclient import AuthClient
from ..client import AsyncClient
from ..config import settings
from ..interface import BaseInterface
from ..interface import PrintInterface
from ..interface import RichJsonInterface
from ..response import check_response

from typing import Dict

async def user_register(
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

async def user_list(client: AuthClient, **kwargs):
    id_list = None
    users = []
    for _id in id_list:
        users.append(await user_show(client, _id))
    return users
    
async def user_show(client: AuthClient, _id: str, **kwargs):
    res = await client.get(f"{settings.FRACTAL_SERVER}/auth/{_id}")
    data = res.json()

    return PrintInterface(
        retcode=0,
        data=data,
    )

async def user_edit(client: AuthClient, _id: str, payload: Dict, **kwargs):
    payload = None # TODO
    res = await client.patch(
        f"{settings.FRACTAL_SERVER}/auth/users/{_id}",
        payload=payload
    )
    data = res.json()

    return PrintInterface(
        retcode=0,
        data=data,
    )

async def user_delete(client: AuthClient, _id: str, **kwargs):
    res = await client.delete(f"{settings.FRACTAL_SERVER}/auth/users/{_id}")
    data = res.json()

    return PrintInterface(
        retcode=0,
        data=data,
    )

async def user_whoami(
    client: AuthClient, **kwargs
) -> PrintInterface:
    res = await client.get(f"{settings.FRACTAL_SERVER}/auth/users/me")
    data = res.json()

    return PrintInterface(
        retcode=0,
        data=data,
    )