from os import environ

import pytest
from devtools import debug

EMAIL_USER = "test@testmail.com"
PWD_USER = "testpassword"


def test_register_as_user(invoke, register_user, caplog):
    with pytest.raises(SystemExit):
        invoke(f"user register {EMAIL_USER} {PWD_USER}")
    debug(caplog.text)
    assert "403" in caplog.text


@pytest.mark.parametrize("is_superuser", [True, False])
def test_register_as_superuser(invoke_as_superuser, is_superuser: bool):
    if is_superuser:
        res = invoke_as_superuser(
            f"user register {EMAIL_USER} {PWD_USER} --superuser"
        )
        debug(res.data)
        assert res.retcode == 0
        assert res.data["is_superuser"]
    else:
        res = invoke_as_superuser(
            f"user register {EMAIL_USER} {PWD_USER} "
            "--slurm-user SOMETHING --cache-dir /absolute --username X"
        )
        debug(res.data)
        assert res.retcode == 0
        assert not res.data["is_superuser"]
    assert res.data["email"] == EMAIL_USER

    # Test that new user is verified (note: for the moment we don't expose the
    # possibility of registering a non-verified user)
    assert res.data["is_verified"]


def test_register_as_superuser_with_batch(invoke_as_superuser):
    # Register a user with the --batch flag
    res = invoke_as_superuser(f"--batch user register {EMAIL_USER} {PWD_USER}")
    user_id = res.data
    debug(user_id)
    assert res.retcode == 0
    # Check that the user exists
    res = invoke_as_superuser(f"user show {user_id}")
    debug(res.data)
    assert res.data["email"] == EMAIL_USER
    assert res.retcode == 0


def test_list_as_user(invoke, register_user, caplog):
    with pytest.raises(SystemExit):
        invoke("user list")
    debug(caplog.text)
    assert "403" in caplog.text


def test_list_as_superuser(invoke_as_superuser, register_user):
    res = invoke_as_superuser("user list")
    debug(res.data)
    assert res.retcode == 0
    list_emails = [user["email"] for user in res.data]
    debug(list_emails)
    assert "admin@fractal.xy" in list_emails
    assert environ["FRACTAL_USER"] in list_emails


def test_show_as_user(invoke, invoke_as_superuser, register_user, caplog):
    # Register a new user
    res = invoke_as_superuser(f"user register {EMAIL_USER} {PWD_USER}")
    user_id = res.data["id"]
    # Call fractal user show
    with pytest.raises(SystemExit):
        invoke(f"user show {user_id}")
    debug(caplog.text)
    assert "403" in caplog.text


def test_show_as_superuser(invoke_as_superuser):
    # Register a new user
    res = invoke_as_superuser(f"user register {EMAIL_USER} {PWD_USER}")
    user_id = res.data["id"]
    # Call fractal user show
    invoke_as_superuser(f"user show {user_id}")
    debug(res.data)
    assert res.retcode == 0
    assert res.data["email"] == EMAIL_USER


def test_edit_as_user(invoke, invoke_as_superuser, register_user, caplog):
    # Register a new user
    res = invoke_as_superuser(f"user register {EMAIL_USER} {PWD_USER}")
    user_id = res.data["id"]
    # Call fractal user edit
    with pytest.raises(SystemExit):
        invoke(f"user edit {user_id} --new-email email@something.xy")
    debug(caplog.text)
    assert "403" in caplog.text


@pytest.mark.parametrize("new_is_superuser", [True, False])
@pytest.mark.parametrize("new_is_non_verified", [True, False])
def test_edit_as_superuser(
    invoke_as_superuser, new_is_superuser, new_is_non_verified
):
    # Register a new user
    res = invoke_as_superuser(f"user register {EMAIL_USER} {PWD_USER}")
    assert res.retcode == 0
    user_id = res.data["id"]
    # Call fractal user edit
    NEW_EMAIL = "asd@asd.new"
    NEW_CACHE_DIR = "/tmp/xxx"
    NEW_SLURM_USER = "new_slurm"
    NEW_USERNAME = "new_username"
    cmd = (
        f"user edit {user_id} "
        f"--new-email {NEW_EMAIL} "
        f"--new-password SOMETHING "
        f"--new-slurm-user {NEW_SLURM_USER} "
        f"--new-username {NEW_USERNAME} "
        f"--new-cache-dir {NEW_CACHE_DIR}"
    )
    if new_is_superuser:
        cmd = f"{cmd} --make-superuser"
    if new_is_non_verified:
        cmd = f"{cmd} --remove-verified"
    debug(cmd)
    res = invoke_as_superuser(cmd)
    debug(res.data)
    assert res.retcode == 0
    assert res.data["email"] == NEW_EMAIL
    assert res.data["cache_dir"] == NEW_CACHE_DIR
    assert res.data["slurm_user"] == NEW_SLURM_USER
    assert res.data["username"] == NEW_USERNAME
    assert res.data["is_superuser"] == new_is_superuser
    if new_is_non_verified:
        assert not res.data["is_verified"]

    BAD_CACHE_DIR = "not_absolute"
    with pytest.raises(SystemExit):
        cmd = f"user edit {user_id} --new-cache-dir {BAD_CACHE_DIR}"
        invoke_as_superuser(cmd)

    # If the user was made a superuser, check that we can go back to normal
    # user
    if new_is_superuser:
        cmd = f"user edit {user_id} --remove-superuser"
        debug(cmd)
        res = invoke_as_superuser(cmd)
        debug(res.data)
        assert res.retcode == 0
        assert not res.data["is_superuser"]

    # If the user was made verified, check that we can go back to normal
    # user
    if new_is_non_verified:
        cmd = f"user edit {user_id} --make-verified"
        debug(cmd)
        res = invoke_as_superuser(cmd)
        debug(res.data)
        assert res.retcode == 0
        assert res.data["is_verified"]


def test_edit_arguments(invoke_as_superuser):
    # Test that superuser flags are mutually exclusive
    with pytest.raises(SystemExit):
        cmd = "user edit SOME_USER_ID --make-superuser --remove-superuser"
        invoke_as_superuser(cmd)


@pytest.mark.skip(
    reason="Delete-user endpoint was removed in fractal-server 1.4.0"
)
def test_delete_as_user(invoke, invoke_as_superuser, register_user, caplog):
    # Register a new user
    res = invoke_as_superuser(f"user register {EMAIL_USER} {PWD_USER}")
    user_id = res.data["id"]
    # Call fractal user edit
    with pytest.raises(SystemExit):
        invoke(f"user delete {user_id}")
    debug(caplog.text)
    assert "403" in caplog.text


@pytest.mark.skip(
    reason="Delete-user endpoint was removed in fractal-server 1.4.0"
)
def test_delete_as_superuser(invoke_as_superuser, caplog):
    # Register a new user
    res = invoke_as_superuser(f"user register {EMAIL_USER} {PWD_USER}")
    user_id = res.data["id"]
    # Call fractal user delete
    res = invoke_as_superuser(f"user delete {user_id}")
    assert res.retcode == 0
    # Check that user was not found
    with pytest.raises(SystemExit):
        res = invoke_as_superuser(f"user show {user_id}")
    debug(caplog.text)
    assert "404" in caplog.text
    assert "Not Found" in caplog.text


def test_whoami_as_user(invoke, register_user):
    res = invoke("user whoami")
    assert res.retcode == 0
    debug(res.data)
    assert res.data["email"] == environ["FRACTAL_USER"]
    assert not res.data["is_superuser"]
    user_id = res.data["id"]

    # Test user whoami with --batch flag
    res = invoke("--batch user whoami")
    debug(res.data)
    assert res.data == user_id
    assert res.retcode == 0


def test_whoami_as_superuser(invoke_as_superuser):
    res = invoke_as_superuser("user whoami")
    assert res.retcode == 0
    debug(res.data)
    assert res.data["email"] == "admin@fractal.xy"
    assert res.data["is_superuser"]
