import io
from datetime import datetime
from datetime import timedelta
from pathlib import Path

import jwt
import pytest

from fractal_client.auth._token_utils import _get_validity_seconds
from fractal_client.auth._token_utils import _is_token_valid
from fractal_client.auth._token_utils import read_and_refresh_token

LONG_KEY = "long-secret-very-secret-key-long"


def test_get_validity_seconds():
    past = (datetime.now() - timedelta(seconds=100)).timestamp()
    future = (datetime.now() + timedelta(seconds=100)).timestamp()

    fake_token = "asdasd"
    assert _get_validity_seconds(fake_token) == -1

    token = jwt.encode({"exp": past}, key=LONG_KEY, algorithm="HS256")
    assert _get_validity_seconds(token) < -100

    token = jwt.encode({"exp": future}, key=LONG_KEY, algorithm="HS256")
    assert _get_validity_seconds(token) > 99


def test_is_token_valid():
    past = (datetime.now() - timedelta(seconds=100)).timestamp()
    near_future = (datetime.now() + timedelta(seconds=1)).timestamp()
    far_future = (datetime.now() + timedelta(seconds=100)).timestamp()

    fake_token = "asdasd"
    assert not _is_token_valid(fake_token)

    token = jwt.encode(payload={"exp": past}, key=LONG_KEY, algorithm="HS256")
    assert not _is_token_valid(token)

    token = jwt.encode({"exp": near_future}, key=LONG_KEY, algorithm="HS256")
    assert not _is_token_valid(token)

    token = jwt.encode({"exp": far_future}, key=LONG_KEY, algorithm="HS256")
    assert _is_token_valid(token)


def test_read_and_refresh_token(tmp_path: Path, monkeypatch):
    invalid_token = jwt.encode(
        payload={"exp": (datetime.now() - timedelta(days=1)).timestamp()},
        key=LONG_KEY,
        algorithm="HS256",
    )
    valid_token = jwt.encode(
        payload={"exp": (datetime.now() + timedelta(days=1)).timestamp()},
        key=LONG_KEY,
        algorithm="HS256",
    )

    valid_path = (tmp_path / "valid").as_posix()
    invalid_path_1 = (tmp_path / "invalid1").as_posix()
    invalid_path_2 = (tmp_path / "invalid2").as_posix()
    invalid_path_3 = (tmp_path / "invalid3").as_posix()
    Path(valid_path).write_text(valid_token)
    Path(invalid_path_1).write_text(invalid_token)
    Path(invalid_path_3).write_text("")

    assert read_and_refresh_token(valid_path) == valid_token

    monkeypatch.setattr("sys.stdin", io.StringIO(valid_token))
    assert read_and_refresh_token(invalid_path_1) == valid_token
    assert Path(invalid_path_1).read_text() == valid_token

    monkeypatch.setattr("sys.stdin", io.StringIO(valid_token))
    assert read_and_refresh_token(invalid_path_2) == valid_token
    assert Path(invalid_path_2).read_text() == valid_token

    with pytest.raises(
        SystemExit,
        match="The token provided is invalid",
    ):
        monkeypatch.setattr("sys.stdin", io.StringIO(invalid_token))
        read_and_refresh_token(invalid_path_3)
    assert Path(invalid_path_3).read_text() == ""
