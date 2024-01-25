from typing import Optional
from typing import Union

from ..authclient import AuthClient
from ..config import settings
from ..interface import PrintInterface
from ..interface import RichJsonInterface
from ..response import check_response


def user_register(
    client: AuthClient,
    *,
    new_email: str,
    new_password: Optional[str] = None,
    slurm_user: Optional[str] = None,
    cache_dir: Optional[str] = None,
    username: Optional[str] = None,
    superuser: bool = False,
    verified: bool = True,  # TODO: this is not currently exposed in the CLI
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

    res = client.post(
        f"{settings.FRACTAL_SERVER}/auth/register/", json=new_user
    )
    data = check_response(res, expected_status_code=201)

    if superuser or verified:
        patch_payload = dict(is_superuser=superuser, is_verified=verified)
        user_id = data["id"]
        res = client.patch(
            f"{settings.FRACTAL_SERVER}/auth/users/{user_id}/",
            json=patch_payload,
        )
        data = check_response(res, expected_status_code=200)

    if batch:
        return PrintInterface(retcode=0, data=data["id"])
    else:
        return RichJsonInterface(retcode=0, data=data)


def user_list(client: AuthClient) -> RichJsonInterface:
    res = client.get(f"{settings.FRACTAL_SERVER}/auth/users/")
    users = check_response(res, expected_status_code=200)
    return RichJsonInterface(retcode=0, data=users)


def user_show(client: AuthClient, *, user_id: str) -> RichJsonInterface:
    res = client.get(f"{settings.FRACTAL_SERVER}/auth/users/{user_id}/")
    user = check_response(res, expected_status_code=200)
    return RichJsonInterface(retcode=0, data=user)


def user_edit(
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
    make_verified: bool = False,
    remove_verified: bool = False,
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
    if make_verified:
        user_update["is_verified"] = True
    if remove_verified:
        user_update["is_verified"] = False
    if new_cache_dir is not None:
        user_update["cache_dir"] = new_cache_dir
    if new_slurm_user is not None:
        user_update["slurm_user"] = new_slurm_user
    if new_username is not None:
        user_update["username"] = new_username

    if not user_update:
        return PrintInterface(retcode=1, data="Nothing to update")
    res = client.patch(
        f"{settings.FRACTAL_SERVER}/auth/users/{user_id}/", json=user_update
    )
    new_user = check_response(res, expected_status_code=200)

    return RichJsonInterface(retcode=0, data=new_user)


def user_whoami(
    client: AuthClient, *, batch: bool
) -> Union[RichJsonInterface, PrintInterface]:
    res = client.get(f"{settings.FRACTAL_SERVER}/auth/current-user/")
    user = check_response(res, expected_status_code=200)
    if batch:
        return PrintInterface(retcode=0, data=user["id"])
    else:
        return RichJsonInterface(retcode=0, data=user)
