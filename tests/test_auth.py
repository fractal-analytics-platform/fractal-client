import pytest
from devtools import debug

from .fixtures_testserver import TESTER
from fractal_client.authclient import AuthenticationError
from fractal_client.authclient import AuthToken


def test_auth_registered(client):
    """
    GIVEN an existing user
    WHEN fetching a token
    THEN authentication goes through
    """
    auth = AuthToken(
        client, username=TESTER["email"], password=TESTER["password"]
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
