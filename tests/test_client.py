from os import environ

import httpx
from devtools import debug

from fractal.client.config import __VERSION__
from fractal.client.newclient import main


DEFAULT_TEST_EMAIL = environ["FRACTAL_USER"]


async def test_version(clisplit, testserver):
    iface = await main(clisplit("version"))
    debug(iface.output)
    assert f"version: {__VERSION__}" in iface.output
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
