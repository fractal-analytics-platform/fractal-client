from os import environ

import pytest
from devtools import debug

# user register


async def test_register(invoke, register_user, caplog):
    EMAIL = "test@testmail.com"
    PASSWORD = "testpassword"
    with pytest.raises(SystemExit):
        await invoke(f"user register {EMAIL} {PASSWORD}")
    debug(caplog.text)
    assert "403" in caplog.text
    assert "Forbidden" in caplog.text


@pytest.mark.parametrize("is_superuser", [True, False])
async def test_register_superuser(invoke_as_superuser, is_superuser: bool):
    EMAIL = "test@testmail.com"
    PASSWORD = "testpassword"
    if is_superuser:
        res = await invoke_as_superuser(
            f"user register {EMAIL} {PASSWORD} --is-superuser"
        )
        debug(res.data)
        assert res.retcode == 0
        assert res.data["is_superuser"]
    else:
        res = await invoke_as_superuser(f"user register {EMAIL} {PASSWORD}")
        debug(res.data)
        assert res.retcode == 0
        assert not res.data["is_superuser"]
    assert res.data["email"] == EMAIL


# user list


async def test_list(invoke, register_user, caplog):
    with pytest.raises(SystemExit):
        await invoke("user list")
    debug(caplog.text)
    assert "403" in caplog.text
    assert "Forbidden" in caplog.text


async def test_list_superuser(invoke_as_superuser):
    res = await invoke_as_superuser("user list")
    assert res.retcode == 0
    debug(res.data)


# user show


async def test_show(invoke, invoke_as_superuser, register_user, caplog):

    # Register a new user
    EMAIL = "test@testmail.com"
    PASSWORD = "testpassword"
    res = await invoke_as_superuser(f"user register {EMAIL} {PASSWORD}")
    user_id = res.data["id"]

    # Call fractal user show
    with pytest.raises(SystemExit):
        await invoke(f"user show {user_id}")
    debug(caplog.text)
    assert "403" in caplog.text
    assert "Forbidden" in caplog.text


async def test_show_superuser(invoke_as_superuser):

    # Register a new user
    EMAIL = "test@testmail.com"
    PASSWORD = "testpassword"
    res = await invoke_as_superuser(f"user register {EMAIL} {PASSWORD}")
    user_id = res.data["id"]

    # Call fractal user show
    await invoke_as_superuser(f"user show {user_id}")
    debug(res.data)
    assert res.retcode == 0


# user edit


async def test_edit(invoke, register_user):
    # TODO
    pass


async def test_edit_superuser(invoke_as_superuser):
    # TODO
    pass


# user delete


async def test_delete(invoke, register_user):
    # TODO
    pass


async def test_delete_superuser(invoke_as_superuser):
    # TODO
    pass


# user whoami


async def test_whoami(invoke, register_user):
    res = await invoke("user whoami")
    assert res.retcode == 0
    debug(res.data)
    assert res.data["email"] == environ["FRACTAL_USER"]
    assert not res.data["is_superuser"]


async def test_whoami_superuser(invoke_as_superuser):
    res = await invoke_as_superuser("user whoami")
    assert res.retcode == 0
    debug(res.data)
    assert res.data["email"] == "admin@fractal.xy"
    assert res.data["is_superuser"]
