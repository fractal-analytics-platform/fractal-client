import argparse as ap
from datetime import datetime
from datetime import timedelta

import jwt
import pytest

from fractal_client.auth._client import _get_validity_seconds
from fractal_client.auth._client import _is_token_valid
from fractal_client.auth._info import AuthInfo
from fractal_client.auth._info import get_auth_info
from fractal_client.auth.validation import get_cmd_handler
from fractal_client.auth.validation import get_fractal_server

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

    token = jwt.encode({"exp": past}, key=LONG_KEY, algorithm="HS256")
    assert not _is_token_valid(token)

    token = jwt.encode({"exp": near_future}, key=LONG_KEY, algorithm="HS256")
    assert not _is_token_valid(token)

    token = jwt.encode({"exp": far_future}, key=LONG_KEY, algorithm="HS256")
    assert _is_token_valid(token)
