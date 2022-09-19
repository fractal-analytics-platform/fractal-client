from os import environ

import httpx
from devtools import debug

from fractal import __VERSION__


DEFAULT_TEST_EMAIL = environ["FRACTAL_USER"]


async def test_version(invoke, testserver):
    iface = await invoke("version")
    debug(iface.data)
    assert f"version: {__VERSION__}" in iface.data
    assert iface.retcode == 0


async def test_server(testserver):
    """
    GIVEN a testserver
    WHEN it gets called
    THEN it replies
    """
    res = httpx.get("http://localhost:10080/api/alive/")
    debug(res.json())
    assert res.status_code == 200


async def test_user(clear_db, register_user):
    debug(register_user)
    assert register_user["email"] == DEFAULT_TEST_EMAIL
