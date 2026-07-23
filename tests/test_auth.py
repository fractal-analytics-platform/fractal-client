import argparse as ap
import io
from datetime import datetime
from datetime import timedelta
from pathlib import Path

import jwt
import pytest
from httpx2 import Response

from fractal_client.auth._client import AuthClient
from fractal_client.auth._client import AuthenticationError
from fractal_client.auth._client import _get_validity_seconds
from fractal_client.auth._client import _is_token_valid
from fractal_client.auth._client import _read_and_refresh_token
from fractal_client.auth._info import AuthInfo
from fractal_client.auth._info import get_auth_info
from fractal_client.auth.validation import get_cmd_handler
from fractal_client.auth.validation import get_fractal_server
from fractal_client.config import settings

LONG_KEY = "long-secret-very-secret-key-long"


def test_get_fractal_server(monkeypatch):
    args = ap.Namespace(fractal_server=None)

    # Success due to env variable
    assert get_fractal_server(args) is not None

    import fractal_client.main

    monkeypatch.setattr(
        fractal_client.main.settings,
        "FRACTAL_SERVER",
        None,
    )
    with pytest.raises(
        SystemExit,
        match="fractal-server URL",
    ):
        get_fractal_server(args)


def test_get_cmd_handler(monkeypatch):
    with pytest.raises(
        SystemExit,
        match="Help message",
    ):
        get_cmd_handler(
            args=ap.Namespace(cmd=None),
            parser_help="Help message",
        )

    with pytest.raises(
        SystemExit,
        match="internal error",
    ):
        get_cmd_handler(
            args=ap.Namespace(cmd="missing-command"),
            parser_help="Help message",
        )

    cmd_handler = get_cmd_handler(
        args=ap.Namespace(cmd="version"),
        parser_help="Help message",
    )
    assert cmd_handler is not None


def test_AuthInfo():
    assert AuthInfo(
        user="a",
        password="a",
        token_path=None,
    ).use_basic_auth
    assert not AuthInfo(
        user=None,
        password=None,
        token_path="/a",
    ).use_basic_auth


def test_get_auth_info():
    with pytest.raises(
        SystemExit,
        match="Invalid authentication",
    ):
        get_auth_info(ap.Namespace(user=None, password=None, token_path=None))

    with pytest.raises(
        SystemExit,
        match="Invalid authentication",
    ):
        get_auth_info(ap.Namespace(user="user", password=None, token_path=None))

    assert get_auth_info(
        ap.Namespace(user=None, password="password", token_path=None)
    ).use_basic_auth

    with pytest.raises(
        SystemExit,
        match="Invalid authentication",
    ):
        get_auth_info(
            ap.Namespace(
                user="user",
                password="password",
                token_path="/token",
            )
        )


def test_get_validity_seconds():
    past = (datetime.now() - timedelta(seconds=100)).timestamp()
    future = (datetime.now() + timedelta(seconds=100)).timestamp()

    fake_token = "asdasd"
    assert _get_validity_seconds(fake_token) == -1

    token = jwt.encode({"exp": past}, key=LONG_KEY, algorithm="HS256")
    assert _get_validity_seconds(token) < -100

    token = jwt.encode({"exp": future}, key=LONG_KEY, algorithm="HS256")
    assert _get_validity_seconds(token) > 99


def test_is_token_valid():
    past = (datetime.now() - timedelta(seconds=100)).timestamp()
    near_future = (datetime.now() + timedelta(seconds=1)).timestamp()
    far_future = (datetime.now() + timedelta(seconds=100)).timestamp()

    fake_token = "asdasd"
    assert not _is_token_valid(fake_token)

    token = jwt.encode(payload={"exp": past}, key=LONG_KEY, algorithm="HS256")
    assert not _is_token_valid(token)

    token = jwt.encode({"exp": near_future}, key=LONG_KEY, algorithm="HS256")
    assert not _is_token_valid(token)

    token = jwt.encode({"exp": far_future}, key=LONG_KEY, algorithm="HS256")
    assert _is_token_valid(token)


def test_read_and_refresh_token(tmp_path: Path, monkeypatch):
    invalid_token = jwt.encode(
        payload={"exp": (datetime.now() - timedelta(days=1)).timestamp()},
        key=LONG_KEY,
        algorithm="HS256",
    )
    valid_token = jwt.encode(
        payload={"exp": (datetime.now() + timedelta(days=1)).timestamp()},
        key=LONG_KEY,
        algorithm="HS256",
    )

    valid_path = (tmp_path / "valid").as_posix()
    invalid_path_1 = (tmp_path / "invalid1").as_posix()
    invalid_path_2 = (tmp_path / "invalid2").as_posix()
    invalid_path_3 = (tmp_path / "invalid3").as_posix()
    Path(valid_path).write_text(valid_token)
    Path(invalid_path_1).write_text(invalid_token)
    Path(invalid_path_3).write_text("")

    assert _read_and_refresh_token(valid_path) == valid_token

    monkeypatch.setattr("sys.stdin", io.StringIO(valid_token))
    assert _read_and_refresh_token(invalid_path_1) == valid_token
    assert Path(invalid_path_1).read_text() == valid_token

    monkeypatch.setattr("sys.stdin", io.StringIO(valid_token))
    assert _read_and_refresh_token(invalid_path_2) == valid_token
    assert Path(invalid_path_2).read_text() == valid_token

    with pytest.raises(
        SystemExit,
        match="The token provided is invalid",
    ):
        monkeypatch.setattr("sys.stdin", io.StringIO(invalid_token))
        _read_and_refresh_token(invalid_path_3)
    assert Path(invalid_path_3).read_text() == ""


def test_response_to_token():
    auth_client = AuthClient(
        fractal_server="http://localhost:8765",
        auth_info=AuthInfo(
            user="fake",
            password="fake",
            token_path=None,
        ),
    )

    with pytest.raises(
        AuthenticationError,
        match="JSONDecodeError",
    ):
        auth_client._response_to_token(Response(status_code=200, content="non-json"))

    with pytest.raises(
        AuthenticationError,
        match="KeyError",
    ):
        auth_client._response_to_token(Response(status_code=200, json={}))


def test_AuthClient(tester, tmp_path: Path, override_settings):

    # Make sure that FRACTAL_CACHE_PATH is in a temporary folder
    override_settings()

    with AuthClient(
        fractal_server="http://localhost:8765",
        auth_info=AuthInfo(
            user=tester["email"],
            password=tester["password"],
            token_path=None,
        ),
    ) as client:
        assert _is_token_valid(client.token)
        valid_token = client.token

    with pytest.raises(AuthenticationError):
        with AuthClient(
            fractal_server="http://localhost:8765",
            auth_info=AuthInfo(
                user=tester["email"],
                password="wrong-password",
                token_path=None,
            ),
        ) as client:
            pass

    custom_token_path = (tmp_path / "token").as_posix()
    Path(custom_token_path).write_text(valid_token)
    with AuthClient(
        fractal_server="http://localhost:8765",
        auth_info=AuthInfo(
            user=None,
            password=None,
            token_path=custom_token_path,
        ),
    ) as client:
        assert _is_token_valid(client.token)
        assert client.token == valid_token

    Path(settings.default_token_path).write_text(valid_token)
    with AuthClient(
        fractal_server="http://localhost:8765",
        auth_info=AuthInfo(
            user=None,
            password=None,
            token_path=None,
        ),
    ) as client:
        assert _is_token_valid(client.token)
        assert client.token == valid_token
