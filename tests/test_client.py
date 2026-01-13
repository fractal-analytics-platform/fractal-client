import shlex
from pathlib import Path

import httpx
import pytest
from devtools import debug
from fractal_client import __VERSION__
from fractal_client.authclient import AuthClient
from fractal_client.client import _verify_authentication_branch
from fractal_client.client import handle
from fractal_client.cmd import version


def test_debug(invoke):
    res = invoke("--debug version")
    assert res.retcode == 0
    debug(res.data)


def test_version(invoke):
    iface = invoke("version")
    debug(iface.data)
    assert f"version: {__VERSION__}" in iface.data
    assert iface.retcode == 0


def test_version_connect_error():
    iface = version("http://localhost:9999")
    debug(iface.data)
    assert f"version: {__VERSION__}" in iface.data
    assert "refused" in iface.data
    assert iface.retcode == 0


def test_server_is_up():
    """
    GIVEN a testserver
    WHEN it gets called
    THEN it replies
    """
    res = httpx.get("http://localhost:8765/api/alive/")
    debug(res.json())
    assert res.status_code == 200


def test_user_whoami(tester, invoke):
    res = invoke("user whoami")
    user = res.data
    debug(user)
    assert res.retcode == 0
    assert user["email"] == tester["email"]


def test_user_override(user_factory, invoke):
    """
    GIVEN a user whose credentials differ from those of the environment
    WHEN the client is invoked with -u and -p
    THEN the credentials are overridden
    """
    EMAIL = "other_user@exact-lab.it"
    PASSWORD = "other_password"
    user_factory(email=EMAIL, password=PASSWORD)

    res = invoke(f"-u {EMAIL} -p {PASSWORD} project list")
    assert res.retcode == 0


def test_bad_credentials(invoke):
    """
    GIVEN a registered user
    WHEN wrong credentials are passed
    THEN the client returns an error
    """
    res = invoke("-u nouser@exact-lab.it -p nopassword project list")
    res.show()
    assert res.retcode != 0
    assert "BAD_CREDENTIALS" in res.data


def test_connecterror(override_settings):
    override_settings(
        FRACTAL_USER="admin@example.org",
        FRACTAL_PASSWORD="1234",
        FRACTAL_SERVER="http://localhost:12345",
    )
    res = handle(shlex.split("fractal user whoami"))
    debug(res.data)
    assert "ConnectError" in res.data
    assert "Hint: is http://localhost:12345 alive?" in res.data


def test_argparse_abbreviation(invoke_as_superuser):
    """
    Check that argparse abbreviations are disabled on at least one command.

    Refs:
    * https://github.com/fractal-analytics-platform/fractal/issues/440
    * https://docs.python.org/3/library/argparse.html#prefix-matching
    """

    # Successful invoke
    res = invoke_as_superuser(
        "user register test@mail.com secret /project-dir --superuser"
    )
    res.show()
    assert res.retcode == 0

    # Failed (abbreviation-based) invoke
    with pytest.raises(SystemExit):
        invoke_as_superuser(
            "user register test2@mail.com secret2 /project-dir --super"
        )


def test_unit_verify_authentication_branch():
    # Valid cases
    _verify_authentication_branch(
        username="xxx",
        password="xxx",
        token_path=None,
    )
    _verify_authentication_branch(
        username=None,
        password=None,
        token_path="xxx",
    )

    # Invalid cases
    for username, password, token_path in [
        (None, None, None),
        ("xx", None, None),
        (None, "xx", None),
        ("xx", "xx", "xx"),
        ("xx", None, "xx"),
        (None, "xx", "xx"),
    ]:
        with pytest.raises(
            ValueError,
            match="Invalid authentication credentials",
        ):
            _verify_authentication_branch(
                username=username,
                password=password,
                token_path=token_path,
            )


def test_invalid_credentials(monkeypatch):
    import fractal_client.client

    monkeypatch.setattr(
        fractal_client.client.settings, "FRACTAL_USER", "some-user"
    )
    monkeypatch.setattr(
        fractal_client.client.settings, "FRACTAL_PASSWORD", None
    )
    interface = handle(shlex.split("fractal user whoami"))
    assert "Invalid authentication credentials" in interface.data
    assert interface.retcode == 1


def test_invalid_token_path():
    cmd = "fractal --token-path missingfile user whoami"
    interface = handle(shlex.split(cmd))
    interface.show()
    assert interface.retcode == 1


def test_valid_token_path(
    tmp_path: Path,
    monkeypatch,
    tester,
):
    # Get valid token
    with AuthClient(
        fractal_server="http://localhost:8765",
        username=tester["email"],
        password=tester["password"],
        token=None,
    ) as client:
        token_data = client.token
        debug(token_data)
    token_path = (tmp_path / "token").as_posix()

    import fractal_client.client

    monkeypatch.setattr(
        fractal_client.client.settings,
        "FRACTAL_SERVER",
        "http://localhost:8765",
    )

    # Use valid token
    with open(token_path, "w") as f:
        f.write(token_data)
    cmd = f"fractal --token-path {token_path} user whoami"
    interface = handle(shlex.split(cmd))
    assert interface.data["email"] == tester["email"]
    assert interface.retcode == 0

    # Use valid token, with newlines
    with open(token_path, "w") as f:
        f.write(f"\n\n{token_data}\n\n\n")
    cmd = f"fractal --token-path {token_path} user whoami"
    interface = handle(shlex.split(cmd))
    assert interface.data["email"] == tester["email"]
    assert interface.retcode == 0


def test_missing_fractal_server(monkeypatch):
    import fractal_client.client

    monkeypatch.setattr(
        fractal_client.client.settings,
        "FRACTAL_SERVER",
        None,
    )
    interface = handle(shlex.split("fractal user whoami"))
    assert "You should set the fractal-server URL" in interface.data
    assert interface.retcode == 1
