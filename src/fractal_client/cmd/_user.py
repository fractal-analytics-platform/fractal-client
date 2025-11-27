from ..authclient import AuthClient
from ..interface import Interface
from ..response import check_response


def user_register(
    client: AuthClient,
    *,
    new_email: str,
    new_password: str,
    new_project_dir: str,
    superuser: bool = False,
    verified: bool = True,  # TODO: this is not currently exposed in the CLI
    batch: bool = False,
) -> Interface:
    new_user = dict(
        email=new_email,
        password=new_password,
        project_dirs=[new_project_dir],
    )

    res = client.post("auth/register/", json=new_user)
    user_data = check_response(res, expected_status_code=201)

    if superuser or verified:
        patch_payload = dict(is_superuser=superuser, is_verified=verified)
        user_id = user_data["id"]
        res = client.patch(
            f"auth/users/{user_id}/",
            json=patch_payload,
        )
        user_data = check_response(res, expected_status_code=200)

    user_id = user_data["id"]

    if batch:
        return Interface(retcode=0, data=user_data["id"])
    else:
        return Interface(retcode=0, data=user_data)


def user_list(client: AuthClient) -> Interface:
    res = client.get("auth/users/")
    users = check_response(res, expected_status_code=200)
    return Interface(retcode=0, data=users)


def user_show(client: AuthClient, *, user_id: str) -> Interface:
    res = client.get(f"auth/users/{user_id}/")
    user = check_response(res, expected_status_code=200)
    return Interface(retcode=0, data=user)


def user_edit(
    client: AuthClient,
    *,
    user_id: str,
    new_email: str | None = None,
    new_password: str | None = None,
    new_project_dir: str | None = None,
    new_profile_id: str | None = None,
    make_superuser: bool = False,
    remove_superuser: bool = False,
    make_verified: bool = False,
    remove_verified: bool = False,
) -> Interface:
    user_update = dict()
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
    if new_project_dir is not None:
        user_update["project_dir"] = new_project_dir
    if new_profile_id is not None:
        user_update["profile_id"] = new_profile_id
    if make_superuser:
        user_update["is_superuser"] = True
    if remove_superuser:
        user_update["is_superuser"] = False
    if make_verified:
        user_update["is_verified"] = True
    if remove_verified:
        user_update["is_verified"] = False

    res = client.patch(f"auth/users/{user_id}/", json=user_update)
    new_user = check_response(res, expected_status_code=200)

    if new_email is not None:
        # Since `fastapi-users` sets `is_verified` to `False` each time the
        # email is updated, we set `is_verified` as specified by the user.
        res = client.patch(
            f"auth/users/{user_id}/",
            json=dict(is_verified=user_update["is_verified"]),
        )
        new_user = check_response(res, expected_status_code=200)

    return Interface(retcode=0, data=new_user)


def user_set_groups(
    client: AuthClient, *, user_id: int, group_ids: list[int]
) -> Interface:
    res = client.post(
        f"auth/users/{user_id}/set-groups/",
        json=dict(group_ids=group_ids),
    )
    user = check_response(res, expected_status_code=200)
    return Interface(retcode=0, data=user)


def user_whoami(
    client: AuthClient, *, batch: bool, viewer_paths: bool = False
) -> Interface:
    res = client.get("auth/current-user/")
    user = check_response(res, expected_status_code=200)

    if batch:
        return Interface(retcode=0, data=user["id"])

    if viewer_paths:
        res = client.get("auth/current-user/allowed-viewer-paths/")
        returned_viewer_paths = check_response(res, expected_status_code=200)
        return Interface(
            retcode=0,
            data=dict(
                **user,
                viewer_paths=returned_viewer_paths,
            ),
        )
    else:
        return Interface(retcode=0, data=user)
