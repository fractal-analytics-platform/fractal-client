from typing import Dict

from ..authclient import AuthClient
from ..client import AsyncClient
# from ..common.schemas import UserRead  # TODO create schema in common
# from ..common.schemas import UserUpdate  # TODO create schema in common
from ..config import settings
from ..interface import PrintInterface
from ..interface import RichJsonInterface
from ..response import check_response


async def user_register(
    client: AsyncClient,
    email: str,
    slurm_user: str,
    password: str = None,
    superuser: bool = False,
    **kwargs,
) -> RichJsonInterface:
    from getpass import getpass

    if password is None:
        password = getpass()
        confirm_password = getpass("Confirm password: ")
        if password == confirm_password:
            res = await client.post(
                f"{settings.FRACTAL_SERVER}/auth/register",
                json=dict(
                    email=email,
                    slurm_user=slurm_user,
                    password=password,
                    is_superuser=superuser
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


async def user_list(client: AuthClient, **kwargs) -> RichJsonInterface:
    id_list = await client.get(f"{settings.FRACTAL_SERVER}/auth/userlist/")
    users = []
    for _id in id_list:
        users.append((await user_show(client, _id)).data)
    return RichJsonInterface(
        retcode=0,
        data=users,
    )


async def user_show(
    client: AuthClient, user_id: str, **kwargs
) -> RichJsonInterface:
    res = await client.get(f"{settings.FRACTAL_SERVER}/auth/{user_id}")
    user = check_response(
        res,
        expected_status_code=200,
        #coerce=UserRead
    )
    return RichJsonInterface(
        retcode=0,
        data=user.dict(),
    )


async def user_edit(
    client: AuthClient, user_id: str, payload: Dict, **user_update_dict
) -> RichJsonInterface:
    user_update = UserUpdate(**user_update_dict)
    payload = user_update.dict(exclude_unset=True)
    if not payload:
        return PrintInterface(retcode=1, data="Nothing to update")

    res = await client.patch(
        f"{settings.FRACTAL_SERVER}/auth/users/{user_id}", json=payload
    )
    new_user = check_response(
        res,
        expected_status_code=200,
        #coerce=UserRead
    )

    return PrintInterface(
        retcode=0,
        data=new_user.dict(),
    )


async def user_delete(
    client: AuthClient, user_id: str, **kwargs
) -> RichJsonInterface:

    res = await client.delete(
        f"{settings.FRACTAL_SERVER}/auth/users/{user_id}"
    )
    check_response(res, expected_status_code=204)
    return RichJsonInterface(retcode=0, data="")


async def user_whoami(client: AuthClient, **kwargs) -> RichJsonInterface:
    res = await client.get(f"{settings.FRACTAL_SERVER}/auth/users/me")
    user = check_response(res, expected_status_code=200)

    return RichJsonInterface(
        retcode=0,
        data=user,
    )