from typing import Dict
from typing import Optional

from ..authclient import AuthClient
from ..config import settings
from ..interface import PrintInterface
from ..interface import RichJsonInterface
from ..response import check_response

# from ..common.schemas import UserRead  # TODO create schema in common
# from ..common.schemas import UserUpdate  # TODO create schema in common


async def user_register(
    client: AuthClient,
    new_email: str,
    new_password: str = None,
    slurm_user: Optional[str] = None,
    is_superuser: bool = False,
    **kwargs,
) -> RichJsonInterface:
    from getpass import getpass

    if new_password is None:
        new_password = getpass()
        confirm_password = getpass("Confirm password: ")
        if new_password == confirm_password:
            res = await client.post(
                f"{settings.FRACTAL_SERVER}/auth/register",
                json=dict(
                    email=new_email,
                    password=new_password,
                    slurm_user=slurm_user,
                ),
            )

            data = check_response(res, expected_status_code=201)
            iface = RichJsonInterface(retcode=0, data=data)
        else:
            iface = PrintInterface(retcode=1, data="Passwords do not match.")
        return iface

    # FIXME: remove slurm_user
    res = await client.post(
        f"{settings.FRACTAL_SERVER}/auth/register",
        json=dict(
            email=new_email, password=new_password, slurm_user=slurm_user
        ),
    )
    data = check_response(res, expected_status_code=201)

    if is_superuser:
        user_id = res.json()['id']
        res = await client.patch(
            f"{settings.FRACTAL_SERVER}/auth/users/edit/{user_id}",
            json={'is_superuser': True}
        )

    iface = RichJsonInterface(retcode=0, data=data)

    return iface


async def user_list(client: AuthClient, **kwargs) -> RichJsonInterface:
    res = await client.get(f"{settings.FRACTAL_SERVER}/auth/userlist")
    users = check_response(res, expected_status_code=200)
    # users = [UserRead(**user) for user in users]  #FIXME
    return RichJsonInterface(
        retcode=0,
        data=users,
    )


async def user_show(
    client: AuthClient, user_id: str, **kwargs
) -> RichJsonInterface:
    res = await client.get(f"{settings.FRACTAL_SERVER}/auth/{user_id}")
    # FIXME how to handle wrong user_id?
    user = check_response(
        res,
        expected_status_code=200,
        # coerce=UserRead
    )
    return RichJsonInterface(
        retcode=0,
        data=user.dict(),
    )


async def user_edit(
    client: AuthClient, user_id: str, payload: Dict, **user_update_dict
) -> RichJsonInterface:
    # user_update = UserUpdate(**user_update_dict)
    # payload = user_update.dict(exclude_unset=True)
    # FIXME: add UserUpdate
    payload = user_update_dict  # FIXME
    if not payload:
        return PrintInterface(retcode=1, data="Nothing to update")

    res = await client.patch(
        f"{settings.FRACTAL_SERVER}/auth/users/{user_id}", json=payload
    )
    new_user = check_response(
        res,
        expected_status_code=200,
        # coerce=UserRead
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
    res = await client.get(f"{settings.FRACTAL_SERVER}/auth/whoami")
    user = check_response(res, expected_status_code=200)

    return RichJsonInterface(
        retcode=0,
        data=user,
    )
