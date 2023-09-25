from os import environ

import pytest

from fractal_client.authclient import AuthenticationError
from fractal_client.authclient import AuthToken


async def test_auth_fail(client):
    """
    GIVEN no user registered
    WHEN when fetching a token
    THEN authentication error is raised
    """
    with pytest.raises(AuthenticationError):
        auth = AuthToken(
            client,
            username=environ.get("FRACTAL_USER"),
            password=environ.get("FRACTAL_PASSWORD"),
        )
        await auth()


async def test_auth_registered(client, register_user):
    """
    GIVEN no user registered
    WHEN when fetching a token
    THEN authentication error is raised
    """
    auth = AuthToken(
        client,
        username=environ.get("FRACTAL_USER"),
        password=environ.get("FRACTAL_PASSWORD"),
    )
    token = await auth()
    assert token
