from fractal.client.client import version
from devtools import debug
import httpx


async def test_version(cli):
    from fractal.client.config import __VERSION__

    response = await cli.invoke(version)
    debug(response.output)
    assert __VERSION__ in response.output


async def test_server(cli, testserver):
    """
    GIVEN a testserver
    WHEN it gets called
    THEN it replies
    """
    res = httpx.get("http://localhost:10080/api/alive/")
    debug(res.json())
    assert res.status_code == 200


async def test_user(user_factory):
    assert False
