from typing import Optional
from typing import Union

from ..authclient import AuthClient
from ..common.schemas.user import UserCreate
from ..common.schemas.user import UserRead
from ..common.schemas.user import UserUpdate
from ..config import settings
from ..interface import PrintInterface
from ..interface import RichJsonInterface
from ..response import check_response


async def user_register(
    client: AuthClient,
    new_email: str,
    new_password: Optional[str] = None,
    slurm_user: Optional[str] = None,
    superuser: bool = False,
    **kwargs,
) -> Union[RichJsonInterface, PrintInterface]:

    new_user = UserCreate(
        email=new_email, password=new_password, slurm_user=slurm_user
    )

    from getpass import getpass

    if new_password is None:
        new_password = getpass()
        confirm_new_password = getpass("Confirm password: ")
        if new_password == confirm_new_password:
            res = await client.post(
                f"{settings.FRACTAL_SERVER}/auth/register",
                json=new_user.dict(exclude_unset=True),
            )
            data = check_response(res, expected_status_code=201)
            iface = RichJsonInterface(retcode=0, data=data)
        else:
            iface = PrintInterface(retcode=1, data="Passwords do not match.")
        return iface

    res = await client.post(
        f"{settings.FRACTAL_SERVER}/auth/register",
        json=new_user.dict(exclude_unset=True),
    )
    data = check_response(res, expected_status_code=201)

    if superuser:
        user_id = data["id"]
        res = await client.patch(
            f"{settings.FRACTAL_SERVER}/auth/users/{user_id}",
            json={"is_superuser": True},
        )
        data = check_response(res, expected_status_code=200)

    iface = RichJsonInterface(retcode=0, data=data)

    return iface


async def user_list(client: AuthClient, **kwargs) -> RichJsonInterface:
    res = await client.get(f"{settings.FRACTAL_SERVER}/auth/userlist")
    users = check_response(res, expected_status_code=200)
    users = [UserRead(**user) for user in users]
    return RichJsonInterface(
        retcode=0,
        data=users,
    )


async def user_show(
    client: AuthClient, user_id: str, **kwargs
) -> RichJsonInterface:
    res = await client.get(f"{settings.FRACTAL_SERVER}/auth/users/{user_id}")
    user = check_response(res, expected_status_code=200, coerce=UserRead)
    return RichJsonInterface(
        retcode=0,
        data=user.dict(),
    )


async def user_edit(
    client: AuthClient,
    user_id: str,
    new_email: Optional[str] = None,
    new_slurm_user: Optional[str] = None,
    new_is_superuser: bool = False,
    **kwargs,
) -> Union[RichJsonInterface, PrintInterface]:
    user_update = UserUpdate(
        email=new_email,
        slurm_user=new_slurm_user,
    )
    if new_is_superuser:
        user_update.is_superuser = True
    payload = user_update.dict(exclude_unset=True)
    if not payload:
        return PrintInterface(retcode=1, data="Nothing to update")

    res = await client.patch(
        f"{settings.FRACTAL_SERVER}/auth/users/{user_id}", json=payload
    )
    new_user = check_response(res, expected_status_code=200, coerce=UserRead)

    return RichJsonInterface(
        retcode=0,
        data=new_user.dict(),
    )


async def user_delete(
    client: AuthClient, user_id: str, **kwargs
) -> PrintInterface:

    res = await client.delete(
        f"{settings.FRACTAL_SERVER}/auth/users/{user_id}"
    )
    check_response(res, expected_status_code=204)
    return PrintInterface(retcode=0, data="")


async def user_whoami(client: AuthClient, **kwargs) -> RichJsonInterface:
    res = await client.get(f"{settings.FRACTAL_SERVER}/auth/whoami")
    user = check_response(res, expected_status_code=200)

    return RichJsonInterface(
        retcode=0,
        data=user,
    )
