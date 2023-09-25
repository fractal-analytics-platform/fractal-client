from typing import Optional
from typing import Union

from ..authclient import AuthClient
from ..config import settings
from ..interface import PrintInterface
from ..interface import RichJsonInterface
from ..response import check_response


async def user_register(
    client: AuthClient,
    *,
    new_email: str,
    new_password: Optional[str] = None,
    slurm_user: Optional[str] = None,
    cache_dir: Optional[str] = None,
    username: Optional[str] = None,
    superuser: bool = False,
    batch: bool = False,
) -> Union[RichJsonInterface, PrintInterface]:

    new_user = dict(
        email=new_email,
        password=new_password,
    )
    if slurm_user:
        new_user["slurm_user"] = slurm_user
    if cache_dir:
        new_user["cache_dir"] = cache_dir
    if username:
        new_user["username"] = username

    from getpass import getpass

    if new_password is None:
        new_password = getpass()
        confirm_new_password = getpass("Confirm password: ")
        if new_password == confirm_new_password:
            new_user["password"] = new_password
        else:
            return PrintInterface(retcode=1, data="Passwords do not match.")

    res = await client.post(
        f"{settings.FRACTAL_SERVER}/auth/register", json=new_user
    )
    data = check_response(res, expected_status_code=201)

    if superuser:
        user_id = data["id"]
        res = await client.patch(
            f"{settings.FRACTAL_SERVER}/auth/users/{user_id}",
            json={"is_superuser": True},
        )
        data = check_response(res, expected_status_code=200)

    if batch:
        return PrintInterface(retcode=0, data=data["id"])
    else:
        return RichJsonInterface(retcode=0, data=data)


async def user_list(client: AuthClient) -> RichJsonInterface:
    res = await client.get(f"{settings.FRACTAL_SERVER}/auth/userlist")
    users = check_response(res, expected_status_code=200)
    return RichJsonInterface(retcode=0, data=users)


async def user_show(client: AuthClient, *, user_id: str) -> RichJsonInterface:
    res = await client.get(f"{settings.FRACTAL_SERVER}/auth/users/{user_id}")
    user = check_response(res, expected_status_code=200)
    return RichJsonInterface(retcode=0, data=user)


async def user_edit(
    client: AuthClient,
    *,
    user_id: str,
    new_email: Optional[str] = None,
    new_password: Optional[str] = None,
    new_slurm_user: Optional[str] = None,
    new_cache_dir: Optional[str] = None,
    new_username: Optional[str] = None,
    make_superuser: bool = False,
    remove_superuser: bool = False,
) -> Union[RichJsonInterface, PrintInterface]:

    user_update = dict()
    if new_email is not None:
        user_update["email"] = new_email
    if new_password is not None:
        user_update["password"] = new_password
    if make_superuser:
        user_update["is_superuser"] = True
    if remove_superuser:
        user_update["is_superuser"] = False
    if new_cache_dir is not None:
        user_update["cache_dir"] = new_cache_dir
    if new_slurm_user is not None:
        user_update["slurm_user"] = new_slurm_user
    if new_username is not None:
        user_update["username"] = new_username

    if not user_update:
        return PrintInterface(retcode=1, data="Nothing to update")
    res = await client.patch(
        f"{settings.FRACTAL_SERVER}/auth/users/{user_id}", json=user_update
    )
    from devtools import debug

    debug(user_update)
    debug(res)
    debug(res.json())
    new_user = check_response(res, expected_status_code=200)

    return RichJsonInterface(retcode=0, data=new_user)


async def user_delete(client: AuthClient, *, user_id: str) -> PrintInterface:

    res = await client.delete(
        f"{settings.FRACTAL_SERVER}/auth/users/{user_id}"
    )
    check_response(res, expected_status_code=204)
    return PrintInterface(retcode=0, data="")


async def user_whoami(
    client: AuthClient, *, batch: bool
) -> Union[RichJsonInterface, PrintInterface]:
    res = await client.get(f"{settings.FRACTAL_SERVER}/auth/whoami")
    user = check_response(res, expected_status_code=200)

    if batch:
        return PrintInterface(retcode=0, data=user["id"])
    else:
        return RichJsonInterface(retcode=0, data=user)
