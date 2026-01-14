import multiprocessing
import os
import shlex
import urllib.request
from datetime import datetime
from os import environ
from pathlib import Path

import pytest
from fractal_client.interface import Interface

# This variable must be defined before the first import of config.py
environ["FRACTAL_SERVER"] = "http://127.0.0.1:8765"
from fractal_client.client import handle  # noqa: E402

# set_start_method("fork") necessary to run tests on MacOS
# https://github.com/pytest-dev/pytest-flask/issues/104#issuecomment-577908228
# https://docs.python.org/3/library/multiprocessing.html#multiprocessing.get_start_method
multiprocessing.set_start_method("fork")


@pytest.fixture(scope="session")
def fractal_tasks_mock(tmpdir_factory) -> str:
    url = (
        "https://github.com/fractal-analytics-platform/"
        "fractal-server/raw/main/tests/v2/fractal_tasks_mock/dist/"
        "fractal_tasks_mock-0.0.1-py3-none-any.whl"
    )
    local_path = (
        tmpdir_factory.mktemp("fractal_task_mock")
        / "fractal_tasks_mock-0.0.1-py3-none-any.whl"
    )
    urllib.request.urlretrieve(url, local_path)
    return str(local_path)


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


@pytest.fixture(scope="session")
def tester(new_name):
    return dict(
        email=f"client_tester_{new_name()}@example.org",
        password="pytest",
        project_dir="/tester-project-dir",
    )


@pytest.fixture
def invoke(tester):
    def __invoke(args: str) -> Interface:
        new_args = (
            f"--user {tester['email']} --password {tester['password']} {args}"
        )
        return handle(_clisplit(new_args))

    return __invoke


@pytest.fixture
def invoke_as_superuser():
    def __invoke(args: str) -> Interface:
        new_args = f"--user admin@example.org --password 1234 {args}"
        return handle(_clisplit(new_args))

    return __invoke


@pytest.fixture
def invoke_as_custom_user():
    def __invoke(args: str, email: str, password: str):
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
    class Counter:
        ind: int = 0

        def __next__(self) -> str:
            self.ind = self.ind + 1
            return (
                "name-"
                f"{self.ind - 1}-"
                f"{os.getpid()}"
                f"{int(datetime.now().timestamp())}"
            )

    names = Counter()

    return lambda: next(names)


from .fixtures_testserver import *  # noqa: 401
