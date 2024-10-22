import pytest
from devtools import debug

from fractal_client.authclient import AuthenticationError
from fractal_client.authclient import AuthToken


def test_auth_registered(client, tester):
    """
    GIVEN an existing user
    WHEN fetching a token
    THEN authentication goes through
    """
    auth = AuthToken(
        client, username=tester["email"], password=tester["password"]
    )
    token = auth()
    assert token


def test_auth_fail(client):
    """
    GIVEN no user registered
    WHEN fetching a token
    THEN authentication error is raised
    """
    with pytest.raises(AuthenticationError) as err:
        auth = AuthToken(client, username="foo", password="bar")
        auth()
    debug(err.value)
