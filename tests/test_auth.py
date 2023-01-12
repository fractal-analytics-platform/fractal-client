from os import environ

import pytest
from devtools import debug

from fractal.authclient import AuthenticationError
from fractal.authclient import AuthToken


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
            slurm_user=environ.get("SLURM_USER"),
        )
        await auth()


@pytest.mark.parametrize("slurm_user", [None, "env", "some_other_slurm_user"])
async def test_auth_registered(client, register_user, slurm_user):
    """
    GIVEN a registered user
    WHEN when fetching a token (either with/without providing the slurm_user)
    THEN token is obtained
    """
    if slurm_user == "env":
        slurm_user = environ.get("SLURM_USER")
    debug(slurm_user)
    auth = AuthToken(
        client,
        username=environ.get("FRACTAL_USER"),
        password=environ.get("FRACTAL_PASSWORD"),
        slurm_user=slurm_user,
    )
    token = await auth()
    assert token
    debug(token)
