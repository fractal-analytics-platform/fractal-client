from os import environ

import pytest
from devtools import debug


EMAIL_USER = "test@testmail.com"
PWD_USER = "testpassword"


async def test_register_as_user(invoke, register_user, caplog):
    with pytest.raises(SystemExit):
        await invoke(f"user register {EMAIL_USER} {PWD_USER}")
    debug(caplog.text)
    assert "403" in caplog.text
    assert "Forbidden" in caplog.text


@pytest.mark.parametrize("is_superuser", [True, False])
async def test_register_as_superuser(invoke_as_superuser, is_superuser: bool):
    if is_superuser:
        res = await invoke_as_superuser(
            f"user register {EMAIL_USER} {PWD_USER} --superuser"
        )
        debug(res.data)
        assert res.retcode == 0
        assert res.data["is_superuser"]
    else:
        res = await invoke_as_superuser(
            f"user register {EMAIL_USER} {PWD_USER}"
        )
        debug(res.data)
        assert res.retcode == 0
        assert not res.data["is_superuser"]
    assert res.data["email"] == EMAIL_USER


async def test_list_as_user(invoke, register_user, caplog):
    with pytest.raises(SystemExit):
        await invoke("user list")
    debug(caplog.text)
    assert "403" in caplog.text
    assert "Forbidden" in caplog.text


async def test_list_as_superuser(invoke_as_superuser, register_user):
    res = await invoke_as_superuser("user list")
    debug(res.data)
    assert res.retcode == 0
    list_emails = [user.email for user in res.data]
    debug(list_emails)
    assert "admin@fractal.xy" in list_emails
    assert environ["FRACTAL_USER"] in list_emails


async def test_show_as_user(
    invoke, invoke_as_superuser, register_user, caplog
):
    # Register a new user
    res = await invoke_as_superuser(f"user register {EMAIL_USER} {PWD_USER}")
    user_id = res.data["id"]
    # Call fractal user show
    with pytest.raises(SystemExit):
        await invoke(f"user show {user_id}")
    debug(caplog.text)
    assert "403" in caplog.text
    assert "Forbidden" in caplog.text


async def test_show_as_superuser(invoke_as_superuser):
    # Register a new user
    res = await invoke_as_superuser(f"user register {EMAIL_USER} {PWD_USER}")
    user_id = res.data["id"]
    # Call fractal user show
    await invoke_as_superuser(f"user show {user_id}")
    debug(res.data)
    assert res.retcode == 0
    assert res.data["email"] == EMAIL_USER


async def test_edit_as_user(
    invoke, invoke_as_superuser, register_user, caplog
):
    # Register a new user
    res = await invoke_as_superuser(f"user register {EMAIL_USER} {PWD_USER}")
    user_id = res.data["id"]
    # Call fractal user edit
    with pytest.raises(SystemExit):
        await invoke(f"user edit {user_id} --new-email email@something.xy")
    debug(caplog.text)
    assert "403" in caplog.text
    assert "Forbidden" in caplog.text


@pytest.mark.parametrize("new_is_superuser", [True, False])
async def test_edit_as_superuser(invoke_as_superuser, new_is_superuser):
    # Register a new user
    res = await invoke_as_superuser(f"user register {EMAIL_USER} {PWD_USER}")
    assert res.retcode == 0
    user_id = res.data["id"]
    # Call fractal user edit
    NEW_EMAIL = "asd@asd.new"
    NEW_SLURM_USER = "new_slurm"
    cmd = (
        f"user edit {user_id} "
        f"--new-email {NEW_EMAIL} "
        f"--new-slurm-user {NEW_SLURM_USER} "
    )
    if new_is_superuser:
        cmd = f"{cmd} --make-superuser"
    debug(cmd)
    res = await invoke_as_superuser(cmd)
    debug(res.data)
    assert res.retcode == 0
    assert res.data["email"] == NEW_EMAIL
    assert res.data["slurm_user"] == NEW_SLURM_USER
    assert res.data["is_superuser"] == new_is_superuser

    # If the user was made a superuser, check that we can go back to normal
    # user
    if new_is_superuser:
        cmd = f"user edit {user_id} --remove-superuser"
        debug(cmd)
        res = await invoke_as_superuser(cmd)
        debug(res.data)
        assert res.retcode == 0
        assert not res.data["is_superuser"]


async def test_edit_arguments(invoke_as_superuser):
    # Test that superuser flags are mutually exclusive
    with pytest.raises(SystemExit):
        cmd = "user edit SOME_USER_ID --make-superuser --remove-superuser"
        await invoke_as_superuser(cmd)


async def test_delete_as_user(
    invoke, invoke_as_superuser, register_user, caplog
):
    # Register a new user
    res = await invoke_as_superuser(f"user register {EMAIL_USER} {PWD_USER}")
    user_id = res.data["id"]
    # Call fractal user edit
    with pytest.raises(SystemExit):
        await invoke(f"user delete {user_id}")
    debug(caplog.text)
    assert "403" in caplog.text
    assert "Forbidden" in caplog.text


async def test_delete_as_superuser(invoke_as_superuser, caplog):
    # Register a new user
    res = await invoke_as_superuser(f"user register {EMAIL_USER} {PWD_USER}")
    user_id = res.data["id"]
    # Call fractal user delete
    res = await invoke_as_superuser(f"user delete {user_id}")
    assert res.retcode == 0
    # Check that user was not found
    with pytest.raises(SystemExit):
        res = await invoke_as_superuser(f"user show {user_id}")
    debug(caplog.text)
    assert "404" in caplog.text
    assert "Not Found" in caplog.text


async def test_whoami_as_user(invoke, register_user):
    res = await invoke("user whoami")
    assert res.retcode == 0
    debug(res.data)
    assert res.data["email"] == environ["FRACTAL_USER"]
    assert not res.data["is_superuser"]


async def test_whoami_as_superuser(invoke_as_superuser):
    res = await invoke_as_superuser("user whoami")
    assert res.retcode == 0
    debug(res.data)
    assert res.data["email"] == "admin@fractal.xy"
    assert res.data["is_superuser"]
