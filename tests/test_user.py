from os import environ

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
