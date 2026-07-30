import io
from datetime import datetime
from datetime import timedelta
from pathlib import Path

import jwt
import pytest

from fractal_client.auth._args import AuthInfo
from fractal_client.auth._client import AuthClient
from fractal_client.config import settings


@pytest.fixture(scope="function")
def valid_token(tester) -> str:
    from fractal_client.auth._token_utils import _is_token_valid

    with AuthClient(
        fractal_server="http://localhost:8765",
        auth_info=AuthInfo(
            user=tester["email"],
            password=tester["password"],
            token_path=None,
        ),
    ) as client:
        token = client.token
        assert _is_token_valid(token)
    return token


def test_auth_commands(
    invoke,
    monkeypatch,
    tmp_path: Path,
    valid_token: str,
    override_settings,
    tester: dict,
):
    override_settings()
    custom_token_path = (tmp_path / "custom_token.txt").as_posix()

    # Fully invalid token
    invalid_token = "not-a-token"

    # Partially valid token: it is a valid and not-expiring JWT token, but its signature
    # is based on a key different from the fractal-server one.
    partially_valid_token = jwt.encode(
        payload={"exp": (datetime.now() + timedelta(days=1)).timestamp()},
        key="some-long-jwt-secret-key-which-should-be-long",
        algorithm="HS256",
    )

    for option, token_path in zip(
        (f"--token-path {custom_token_path}", ""),
        (custom_token_path, settings.default_token_path),
    ):
        assert not Path(token_path).exists()

        with pytest.raises(SystemExit, match="not found"):
            invoke(f"{option} auth check-token")

        res = invoke(f"{option} auth clear-token")
        assert res.retcode == 1
        assert "not found" in res.data

        monkeypatch.setattr("sys.stdin", io.StringIO(invalid_token))
        res = invoke(f"{option} auth set-token")
        assert res.retcode == 1
        assert "invalid" in res.data
        assert not Path(token_path).exists()

        monkeypatch.setattr("sys.stdin", io.StringIO(valid_token))
        res = invoke(f"{option} auth set-token")
        assert res.retcode == 0
        assert "Token written" in res.data
        assert Path(token_path).exists()

        res = invoke(f"{option} auth check-token")
        assert f"Valid token for {tester['email']}" in res.data
        assert res.retcode == 0

        res = invoke(f"{option} auth clear-token")
        assert res.retcode == 0
        assert "Removed" in res.data

        monkeypatch.setattr("sys.stdin", io.StringIO(partially_valid_token))
        res = invoke(f"{option} auth set-token")
        assert res.retcode == 0
        assert "Token written" in res.data
        assert Path(token_path).exists()

        res = invoke(f"{option} auth check-token")
        assert res.retcode == 1
        assert "Token verification failed." in res.data
