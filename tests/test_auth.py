import pytest

from fractal.client._auth import AuthenticationError
from fractal.client._auth import AuthToken


async def test_auth_fail(testserver, client):
    """
    GIVEN no user registered
    WHEN when fetching a token
    THEN authentication error is raised
    """
    with pytest.raises(AuthenticationError):
        auth = AuthToken(client)
        await auth()


async def test_auth_registerd(testserver, client, register_user):
    """
    GIVEN no user registered
    WHEN when fetching a token
    THEN authentication error is raised
    """
    auth = AuthToken(client)
    token = await auth()
    assert token
