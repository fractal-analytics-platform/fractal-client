from pathlib import Path

from fractal_client.config import Settings


def test_default_token_path():
    expected_path = (Path.home() / ".cache/fractal" / "token.txt").as_posix()
    assert Settings().default_token_path == expected_path
