import multiprocessing
import shlex
from os import environ
from pathlib import Path

import pytest
from httpx import AsyncClient

TASKS_CACHE_FILENAME = "tasks"


environ["FRACTAL_USER"] = "test@fake-exact-lab.it"
environ["FRACTAL_PASSWORD"] = "password"
environ["FRACTAL_SERVER"] = "http://127.0.0.1:10080"
environ["DB_ECHO"] = "0"

# set_start_method("fork") necessary to run tests on MacOS
# https://github.com/pytest-dev/pytest-flask/issues/104#issuecomment-577908228
# https://docs.python.org/3/library/multiprocessing.html#multiprocessing.get_start_method
multiprocessing.set_start_method("fork")


@pytest.fixture(scope="session")
def testdata_path() -> Path:
    return Path(__file__).parent / "data"


@pytest.fixture(scope="session")
def event_loop():
    import asyncio

    loop = asyncio.get_event_loop()
    yield loop
    loop.close()


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


@pytest.fixture
async def invoke(clisplit):
    from fractal_client.client import handle

    async def __invoke(args: str):
        return await handle(clisplit(args))

    return __invoke


@pytest.fixture
async def invoke_as_superuser(clisplit):
    from fractal_client.client import handle

    async def __invoke(args: str):
        new_args = f"--user admin@fractal.xy --password 1234 {args}"
        return await handle(clisplit(new_args))

    return __invoke


@pytest.fixture
def clear_task_cache():
    from fractal_client.config import settings

    # This is a workaround to clean up the state before the test. The right way
    # would be to inject a new (function-scoped) FRACTAL_CACHE_PATH variable
    # for each test
    cache_dir = Path(settings.FRACTAL_CACHE_PATH)
    cache_file = cache_dir / TASKS_CACHE_FILENAME
    cache_file.unlink(missing_ok=True)
    yield
    cache_file.unlink(missing_ok=True)


from .fixtures_testserver import *  # noqa: 401
