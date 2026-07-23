from pathlib import Path

import pytest
from httpx2 import Response

from fractal_client.auth._args import AuthInfo
from fractal_client.auth._client import AuthClient
from fractal_client.auth._client import AuthenticationError
from fractal_client.auth._token_utils import _is_token_valid
from fractal_client.config import settings


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
        auth_client._get_token_from_response(
            Response(status_code=200, content="non-json")
        )

    with pytest.raises(
        AuthenticationError,
        match="KeyError",
    ):
        auth_client._get_token_from_response(Response(status_code=200, json={}))


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
