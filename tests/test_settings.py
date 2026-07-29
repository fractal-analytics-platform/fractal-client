from pathlib import Path

from fractal_client.config import Settings


def test_default_token_path():
    expected_path = (Path.home() / ".cache/fractal" / "token.txt").as_posix()
    assert Settings().default_token_path == expected_path


def test_fractal_web(monkeypatch):
    clean_url = "https://fractal1.example.org"
    monkeypatch.setenv("FRACTAL_WEB", clean_url)
    assert Settings().FRACTAL_WEB == clean_url

    dirty_url = "https://fractal2.example.org/"
    monkeypatch.setenv("FRACTAL_WEB", dirty_url)
    assert Settings().FRACTAL_WEB == dirty_url[:-1]
