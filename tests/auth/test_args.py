import argparse as ap

import pytest

from fractal_client.auth._args import AuthInfo
from fractal_client.auth._args import get_auth_info
from fractal_client.auth._args import get_fractal_server

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
        get_auth_info(ap.Namespace(user="user", password=None, token_path=None))

    with pytest.raises(
        SystemExit,
        match="Invalid authentication",
    ):
        get_auth_info(ap.Namespace(user=None, password="password", token_path=None))

    assert get_auth_info(
        ap.Namespace(user="user", password="password", token_path=None)
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
