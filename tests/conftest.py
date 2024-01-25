import multiprocessing
import shlex
from os import environ
from pathlib import Path

import pytest
from httpx import Client


# These three variables must be defined before the first import of config.py
environ["FRACTAL_SERVER"] = "http://127.0.0.1:10080"
environ["FRACTAL_USER"] = "test@fake-exact-lab.it"
environ["FRACTAL_PASSWORD"] = "password"

environ["FRACTAL_USERNAME"] = "myusername"


# set_start_method("fork") necessary to run tests on MacOS
# https://github.com/pytest-dev/pytest-flask/issues/104#issuecomment-577908228
# https://docs.python.org/3/library/multiprocessing.html#multiprocessing.get_start_method
multiprocessing.set_start_method("fork")


@pytest.fixture(autouse=True, scope="function")
def clear_cache(tmp_path, monkeypatch):
    """
    Note that this fixture is automatically used in **all** tests.
    """
    import fractal_client.config

    monkeypatch.setattr(
        fractal_client.config.settings,
        "FRACTAL_CACHE_PATH",
        str(tmp_path),
    )


@pytest.fixture(scope="session")
def testdata_path() -> Path:
    return Path(__file__).parent / "data"


@pytest.fixture
def client():
    with Client(timeout=10) as client:
        yield client


@pytest.fixture
def client_superuser():
    from fractal_client.authclient import AuthClient

    with AuthClient(
        username="admin@fractal.xy",
        password="1234",
    ) as client_superuser:
        yield client_superuser


def _clisplit(args: str):
    return shlex.split(f"fractal {args}")


def _remove_session():
    from fractal_client.config import settings

    cache_dir = Path(settings.FRACTAL_CACHE_PATH)
    cache_file = cache_dir / "session"
    cache_file.unlink(missing_ok=True)


@pytest.fixture
def invoke():
    from fractal_client.client import handle

    def __invoke(args: str):
        _remove_session()
        return handle(_clisplit(args))

    return __invoke


@pytest.fixture
def invoke_as_superuser():
    from fractal_client.client import handle

    def __invoke(args: str):
        _remove_session()
        new_args = f"--user admin@fractal.xy --password 1234 {args}"
        return handle(_clisplit(new_args))

    return __invoke


@pytest.fixture(scope="function")
def override_settings(monkeypatch, tmp_path):
    import fractal_client.config

    def _override_settings(
        FRACTAL_CACHE_PATH=str(tmp_path),
        FRACTAL_USER=None,
        FRACTAL_PASSWORD=None,
        FRACTAL_SERVER=None,
    ):
        monkeypatch.setattr(
            fractal_client.config.settings,
            "FRACTAL_CACHE_PATH",
            FRACTAL_CACHE_PATH,
        )
        monkeypatch.setattr(
            fractal_client.config.settings,
            "FRACTAL_USER",
            FRACTAL_USER,
        )
        monkeypatch.setattr(
            fractal_client.config.settings,
            "FRACTAL_PASSWORD",
            FRACTAL_PASSWORD,
        )
        if FRACTAL_SERVER is not None:
            monkeypatch.setattr(
                fractal_client.config.settings,
                "FRACTAL_SERVER",
                FRACTAL_SERVER,
            )

    return _override_settings


from .fixtures_testserver import *  # noqa: 401
