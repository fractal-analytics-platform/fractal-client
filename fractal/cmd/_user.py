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
    cache_dir: Optional[str] = None,
    superuser: bool = False,
    batch: bool = False,
    **kwargs,
) -> Union[RichJsonInterface, PrintInterface]:

    user_dict = dict(
        email=new_email,
        password=new_password,
    )
    if slurm_user:
        user_dict["slurm_user"] = slurm_user
    if cache_dir:
        user_dict["cache_dir"] = cache_dir
    new_user = UserCreate(**user_dict)

    from getpass import getpass

    if new_password is None:
        new_password = getpass()
        confirm_new_password = getpass("Confirm password: ")
        if new_password == confirm_new_password:
            new_user.password = new_password
        else:
            return PrintInterface(retcode=1, data="Passwords do not match.")

    res = await client.post(
        f"{settings.FRACTAL_SERVER}/auth/register",
        json=new_user.dict(exclude_unset=True, exclude_none=True),
    )
    data = check_response(res, expected_status_code=201, coerce=UserRead)

    if superuser:
        user_id = data.id
        res = await client.patch(
            f"{settings.FRACTAL_SERVER}/auth/users/{user_id}",
            json={"is_superuser": True},
        )
        data = check_response(res, expected_status_code=200, coerce=UserRead)

    if batch:
        return PrintInterface(retcode=0, data=data.id)
    else:
        return RichJsonInterface(retcode=0, data=data.dict())


async def user_list(client: AuthClient, **kwargs) -> RichJsonInterface:
    res = await client.get(f"{settings.FRACTAL_SERVER}/auth/userlist")
    users = check_response(res, expected_status_code=200)
    users = [UserRead(**user).dict() for user in users]
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
    new_password: Optional[str] = None,
    new_slurm_user: Optional[str] = None,
    new_cache_dir: Optional[str] = None,
    make_superuser: bool = False,
    remove_superuser: bool = False,
    **kwargs,
) -> Union[RichJsonInterface, PrintInterface]:
    user_dict = dict(
        email=new_email,
        password=new_password,
    )
    if new_cache_dir is not None:
        user_dict["cache_dir"] = new_cache_dir
    if new_slurm_user is not None:
        user_dict["slurm_user"] = new_slurm_user
    user_update = UserUpdate(**user_dict)

    if make_superuser:
        user_update.is_superuser = True
    if remove_superuser:
        user_update.is_superuser = False

    payload = user_update.dict(exclude_unset=True, exclude_none=True)
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


async def user_whoami(
    client: AuthClient, batch: bool, **kwargs
) -> Union[RichJsonInterface, PrintInterface]:
    res = await client.get(f"{settings.FRACTAL_SERVER}/auth/whoami")
    user = check_response(res, expected_status_code=200, coerce=UserRead)

    if batch:
        return PrintInterface(retcode=0, data=user.id)
    else:
        return RichJsonInterface(retcode=0, data=user.dict())
