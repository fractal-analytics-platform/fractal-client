import multiprocessing
import shlex
from os import environ
from pathlib import Path

import pytest
from httpx import AsyncClient


TASKS_CACHE_FILENAME = "tasks"

environ["FRACTAL_USER"] = "test@fake-exact-lab.it"
environ["FRACTAL_PASSWORD"] = "password"
environ["FRACTAL_USERNAME"] = "myusername"
environ["FRACTAL_SERVER"] = "http://127.0.0.1:10080"


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
async def client():
    async with AsyncClient(timeout=10) as client:
        yield client


@pytest.fixture
async def client_superuser():
    from fractal_client.authclient import AuthClient

    async with AuthClient(
        username="admin@fractal.xy",
        password="1234",
    ) as client_superuser:
        yield client_superuser


@pytest.fixture(scope="session")
def clisplit():
    def __clisplit(args: str):
        return shlex.split(f"fractal {args}")

    return __clisplit


def remove_session():
    from fractal_client.config import settings

    cache_dir = Path(settings.FRACTAL_CACHE_PATH)
    cache_file = cache_dir / "session"
    cache_file.unlink(missing_ok=True)


@pytest.fixture
async def invoke(clisplit):
    from fractal_client.client import handle

    async def __invoke(args: str):
        remove_session()

        interface = await handle(clisplit(args))
        return interface

    return __invoke


@pytest.fixture
async def invoke_as_superuser(clisplit):
    from fractal_client.client import handle

    async def __invoke(args: str):
        remove_session()

        new_args = f"--user admin@fractal.xy --password 1234 {args}"
        interface = await handle(clisplit(new_args))
        return interface

    return __invoke


@pytest.fixture(scope="function")
def override_settings(monkeypatch, tmp_path):
    import fractal_client.config

    def _override_settings(
        FRACTAL_CACHE_PATH=str(tmp_path),
        FRACTAL_USER=None,
        FRACTAL_PASSWORD=None,
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

    return _override_settings


from .fixtures_testserver import *  # noqa: 401
