from os import environ

import pytest
from devtools import debug


async def test_whoami(register_user, invoke):
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


async def test_list(register_user, invoke, caplog):
    with pytest.raises(SystemExit):
        await invoke("user list")
    debug(caplog.text)
    assert "403" in caplog.text
    assert "Forbidden" in caplog.text


async def test_list_superuser(invoke_as_superuser):
    res = await invoke_as_superuser("user list")
    assert res.retcode == 0
    debug(res.data)


async def test_show_superuser(invoke_as_superuser):
    res = await invoke_as_superuser("user register a@b.c abc FK")
    debug(res, res.data)
    res = await invoke_as_superuser("user register c@b.a cba AE")
    debug(res, res.data)
