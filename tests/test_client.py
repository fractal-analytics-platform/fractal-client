import shlex
from os import environ

import httpx
import pytest
from devtools import debug

from fractal_client import __VERSION__
from fractal_client.client import MissingCredentialsError


DEFAULT_TEST_EMAIL = environ["FRACTAL_USER"]


async def test_version(invoke):
    iface = await invoke("version")
    debug(iface.data)
    assert f"version: {__VERSION__}" in iface.data
    assert iface.retcode == 0


async def test_server():
    """
    GIVEN a testserver
    WHEN it gets called
    THEN it replies
    """
    res = httpx.get("http://localhost:10080/api/alive/")
    debug(res.json())
    assert res.status_code == 200


async def test_register_user(register_user, invoke):
    res = await invoke("user whoami")
    user = res.data
    debug(user)
    assert res.retcode == 0
    assert user["email"] == DEFAULT_TEST_EMAIL


async def test_user_override(user_factory, invoke):
    """
    GIVEN a user whose credentials differ from those of the environment
    WHEN the client is invoked with -u and -p
    THEN the credentials are overridden
    """
    EMAIL = "other_user@exact-lab.it"
    PASSWORD = "other_password"
    await user_factory(email=EMAIL, password=PASSWORD)

    res = await invoke(f"-u {EMAIL} -p {PASSWORD} project list")
    assert res.retcode == 0


async def test_bad_credentials(invoke):
    """
    GIVEN a registered user
    WHEN wrong credentials are passed
    THEN the client returns an error
    """
    res = await invoke("-u nouser@exact-lab.it -p nopassword project list")
    res.show()
    assert res.retcode != 0
    assert "BAD_CREDENTIALS" in res.data


async def test_missing_credentials(monkeypatch):
    """
    GIVEN an invocation with missing credentials
    THEN the client raises a MissingCredentialsError
    """

    # Define patched settings
    from fractal_client.config import Settings

    patched_settings = Settings()
    patched_settings.FRACTAL_USER = None

    with monkeypatch.context() as m:
        import fractal_client
        from fractal_client.client import handle

        m.setattr(fractal_client.client, "settings", patched_settings)
        debug(fractal_client.config.settings)
        with pytest.raises(MissingCredentialsError) as e:
            await handle(shlex.split("fractal user whoami"))
        debug(e.value)
        debug(e.value.args[0])
        assert "FRACTAL_USER" in e.value.args[0]


async def test_argparse_abbreviation(invoke_as_superuser):
    """
    Check that argparse abbreviations are disabled on at least one command.

    Refs:
    * https://github.com/fractal-analytics-platform/fractal/issues/440
    * https://docs.python.org/3/library/argparse.html#prefix-matching
    """

    # Successful invoke
    res = await invoke_as_superuser(
        "user register test@mail.com secret --superuser"
    )
    res.show()
    assert res.retcode == 0

    # Failed (abbreviation-based) invoke
    with pytest.raises(SystemExit):
        await invoke_as_superuser(
            "user register test2@mail.com secret2 --super"
        )
