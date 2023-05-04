from os import environ
from pathlib import Path

import jwt
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


@pytest.mark.skip("this is not ready")
async def test_wrong_token_in_cache(
    client,
    register_user,
    invoke,
):
    """
    GIVEN an existing cache/session file, with a non-expired token
    WHEN preparing the headers (i.e. calling AuthToken.__call__)
    THEN we reproduce a know error (see issue 380)
    """

    from fractal.config import settings

    # Create invalid token
    cache_dir = Path(f"{settings.FRACTAL_CACHE_PATH}").expanduser()
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_file = cache_dir / "session"
    debug(cache_file)
    wrong_token = jwt.encode({"some": "payload"}, "not_secret")
    debug("REPLACE CACHE WITH ARBITRARY TOKEN")
    with cache_file.open("w") as f:
        f.write(wrong_token)

    # Create AuthToken, and check that it is valid
    auth = AuthToken(
        client,
        username=environ.get("FRACTAL_USER"),
        password=environ.get("FRACTAL_PASSWORD"),
    )
    debug(auth)
    debug(auth.valid)
    assert auth.valid
    token = await auth()
    debug(token)
    assert token

    # Call fractal-server with wrong cached token, this will fail with 401 and
    # then SystemExit
    with pytest.raises(SystemExit) as e:
        await invoke("project new MyProject")
    debug(e.value)

    # FIXME: we are actually trying to reproduce a 500 error, rather than 401

    # FIXME: use tmp_path as FRACTAL_CACHE_PATH, to avoid the following line
    cache_file.unlink()
