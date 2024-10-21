import json
import sys
from pathlib import Path
from typing import Optional

from ..authclient import AuthClient
from ..config import settings
from ..interface import Interface
from ..response import check_response


def user_register(
    client: AuthClient,
    *,
    new_email: str,
    new_password: str,
    slurm_user: Optional[str] = None,
    cache_dir: Optional[str] = None,
    username: Optional[str] = None,
    superuser: bool = False,
    verified: bool = True,  # TODO: this is not currently exposed in the CLI
    batch: bool = False,
) -> Interface:
    new_user = dict(
        email=new_email,
        password=new_password,
    )
    if username:
        new_user["username"] = username

    new_settings = dict()
    if slurm_user:
        new_settings["slurm_user"] = slurm_user
    if cache_dir:
        new_settings["cache_dir"] = cache_dir

    res = client.post(
        f"{settings.FRACTAL_SERVER}/auth/register/", json=new_user
    )
    user_data = check_response(res, expected_status_code=201)

    if superuser or verified:
        patch_payload = dict(is_superuser=superuser, is_verified=verified)
        user_id = user_data["id"]
        res = client.patch(
            f"{settings.FRACTAL_SERVER}/auth/users/{user_id}/",
            json=patch_payload,
        )
        user_data = check_response(res, expected_status_code=200)

    user_id = user_data["id"]
    if new_settings == {}:
        res = client.get(
            f"{settings.FRACTAL_SERVER}/auth/users/{user_id}/settings/"
        )
        user_settings = check_response(res, expected_status_code=200)
    else:
        res = client.patch(
            f"{settings.FRACTAL_SERVER}/auth/users/{user_id}/settings/",
            json=new_settings,
        )
        user_settings = check_response(res, expected_status_code=200)

    if batch:
        return Interface(retcode=0, data=user_data["id"])
    else:
        user_data_with_settings = dict(settings=user_settings, **user_data)
        return Interface(retcode=0, data=user_data_with_settings)


def user_list(client: AuthClient) -> Interface:
    res = client.get(f"{settings.FRACTAL_SERVER}/auth/users/")
    users = check_response(res, expected_status_code=200)
    return Interface(retcode=0, data=users)


def user_show(client: AuthClient, *, user_id: str) -> Interface:
    res = client.get(f"{settings.FRACTAL_SERVER}/auth/users/{user_id}/")
    user = check_response(res, expected_status_code=200)
    user_id = user["id"]
    res = client.get(
        f"{settings.FRACTAL_SERVER}/auth/users/{user_id}/settings/"
    )
    user_settings = check_response(res, expected_status_code=200)
    user_with_settings = dict(settings=user_settings, **user)
    return Interface(retcode=0, data=user_with_settings)


def user_edit(
    client: AuthClient,
    *,
    user_id: str,
    new_email: Optional[str] = None,
    new_password: Optional[str] = None,
    new_username: Optional[str] = None,
    new_slurm_user: Optional[str] = None,
    new_cache_dir: Optional[str] = None,
    new_ssh_settings_json: Optional[str] = None,
    make_superuser: bool = False,
    remove_superuser: bool = False,
    make_verified: bool = False,
    remove_verified: bool = False,
) -> Interface:
    user_update = dict()
    settings_update = dict()
    if new_email is not None:
        if (make_verified is False) and (remove_verified is False):
            # Since `fastapi-users` sets `is_verified` to `False` each time the
            # email is updated, we force the user to make explicit whether the
            # account is verified or not.
            return Interface(
                retcode=1,
                data=(
                    "Cannot use `--new-email` without `--make-verified` or "
                    "`--remove-verified`"
                ),
            )
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
    if new_username is not None:
        user_update["username"] = new_username
    if new_cache_dir is not None:
        settings_update["cache_dir"] = new_cache_dir
    if new_slurm_user is not None:
        settings_update["slurm_user"] = new_slurm_user
    if new_ssh_settings_json is not None:
        new_ssh_settings_json_path = Path(new_ssh_settings_json)
        if not new_ssh_settings_json_path.exists():
            sys.exit(f"Invalid {new_ssh_settings_json=}. File does not exist.")
        with new_ssh_settings_json_path.open("r") as f:
            ssh_settings = json.load(f)
        __ALLOWED_KEYS__ = (
            "ssh_host",
            "ssh_username",
            "ssh_private_key_path",
            "ssh_tasks_dir",
            "ssh_jobs_dir",
        )
        for key, value in ssh_settings.items():
            if key in __ALLOWED_KEYS__:
                settings_update[key] = value
            else:
                sys.exit(f"Invalid {key=} in {new_ssh_settings_json=}.")

    res = client.patch(
        f"{settings.FRACTAL_SERVER}/auth/users/{user_id}/", json=user_update
    )
    new_user = check_response(res, expected_status_code=200)

    if new_email is not None:
        # Since `fastapi-users` sets `is_verified` to `False` each time the
        # email is updated, we set `is_verified` as specified by the user.
        res = client.patch(
            f"{settings.FRACTAL_SERVER}/auth/users/{user_id}/",
            json=dict(is_verified=user_update["is_verified"]),
        )
        new_user = check_response(res, expected_status_code=200)

    if settings_update == {}:
        res = client.get(
            f"{settings.FRACTAL_SERVER}/auth/users/{user_id}/settings/"
        )
        user_settings = check_response(res, expected_status_code=200)
    else:
        res = client.patch(
            f"{settings.FRACTAL_SERVER}/auth/users/{user_id}/settings/",
            json=settings_update,
        )
        user_settings = check_response(res, expected_status_code=200)

    new_user_with_settings = dict(settings=user_settings, **new_user)
    return Interface(retcode=0, data=new_user_with_settings)


def user_whoami(
    client: AuthClient, *, batch: bool, viewer_paths: bool = False
) -> Interface:
    res = client.get(f"{settings.FRACTAL_SERVER}/auth/current-user/")
    user = check_response(res, expected_status_code=200)

    if batch:
        return Interface(retcode=0, data=user["id"])

    res = client.get(f"{settings.FRACTAL_SERVER}/auth/current-user/settings/")
    user_settings = check_response(res, expected_status_code=200)
    user_with_settings = dict(**user, settings=user_settings)

    if viewer_paths:
        res = client.get(
            f"{settings.FRACTAL_SERVER}/auth/current-user/viewer-paths/"
        )
        returned_viewer_paths = check_response(res, expected_status_code=200)
        return Interface(
            retcode=0,
            data=dict(
                **user_with_settings,
                viewer_paths=returned_viewer_paths,
            ),
        )
    else:
        return Interface(retcode=0, data=user_with_settings)
