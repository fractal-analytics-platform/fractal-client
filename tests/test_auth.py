from os import environ

import pytest

from fractal_client.authclient import AuthenticationError
from fractal_client.authclient import AuthToken


def test_auth_registered(client, register_user):
    """
    GIVEN an existing user
    WHEN fetching a token
    THEN authentication goes through
    """
    auth = AuthToken(
        client,
        username=environ.get("FRACTAL_USER"),
        password=environ.get("FRACTAL_PASSWORD"),
    )
    token = auth()
    assert token


def test_auth_fail(client):
    """
    GIVEN no user registered
    WHEN fetching a token
    THEN authentication error is raised
    """
    with pytest.raises(AuthenticationError):
        auth = AuthToken(client, username="foo", password="bar")
        auth()
