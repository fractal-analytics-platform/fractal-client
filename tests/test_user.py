from os import environ

import httpx
from devtools import debug

from fractal import __VERSION__


async def test_whoami(register_user, invoke):
    res = await invoke("user whoami")
    assert res.retcode == 0
    debug(res.data)
    assert res.data["email"] == environ["FRACTAL_USER"]
    assert not res.data["is_superuser"]

async def test_whoami_superuser(default_superuser, invoke):
    res = await invoke("user whoami")
    assert res.retcode == 0
    debug(res.data)
    assert res.data["email"] == "admin@fractal.xy"
    assert res.data["is_superuser"]

    