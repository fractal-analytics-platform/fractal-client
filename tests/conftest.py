import multiprocessing
import shlex
from os import environ
from pathlib import Path

import pytest


# This variable must be defined before the first import of config.py
environ["FRACTAL_SERVER"] = "http://127.0.0.1:8765"
from fractal_client.client import handle  # noqa: E402

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


def _clisplit(args: str):
    return shlex.split(f"fractal {args}")


def _remove_session():
    from fractal_client.config import settings

    cache_dir = Path(settings.FRACTAL_CACHE_PATH)
    cache_file = cache_dir / "session"
    cache_file.unlink(missing_ok=True)


@pytest.fixture(scope="session")
def tester():
    return dict(email="client_tester@example.org", password="pytest")


@pytest.fixture
def invoke(tester):
    def __invoke(args: str):
        _remove_session()
        new_args = (
            f"--user {tester['email']} --password {tester['password']} {args}"
        )
        return handle(_clisplit(new_args))

    return __invoke


@pytest.fixture
def invoke_as_superuser():
    def __invoke(args: str):
        _remove_session()
        new_args = f"--user admin@fractal.xy --password 1234 {args}"
        return handle(_clisplit(new_args))

    return __invoke


@pytest.fixture
def invoke_as_custom_user():
    def __invoke(args: str, email: str, password: str):
        _remove_session()
        new_args = f"--user {email} --password {password} {args}"
        return handle(_clisplit(new_args))

    return __invoke


@pytest.fixture
def superuser(invoke_as_superuser):
    return invoke_as_superuser("user whoami").data


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


@pytest.fixture(scope="session")
def new_name():
    class Counter(object):
        ind: int = 0

        def __next__(self):
            self.ind = self.ind + 1
            return f"name{self.ind - 1}"

    names = Counter()

    return lambda: next(names)


from .fixtures_testserver import *  # noqa: 401
