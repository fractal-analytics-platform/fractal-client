import os
import shlex
from os import environ
from pathlib import Path

import pytest
from httpx import AsyncClient


environ["FRACTAL_USER"] = "test@fake-exact-lab.it"
environ["FRACTAL_PASSWORD"] = "password"
environ["FRACTAL_SERVER"] = "http://127.0.0.1:10080"
environ["DB_ECHO"] = "0"
environ["SLURM_USER"] = "slurm_user"


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
    async with AsyncClient() as client:
        yield client


@pytest.fixture(scope="session")
def clisplit():
    def __clisplit(args: str):
        return shlex.split(f"fractal {args}")

    return __clisplit


@pytest.fixture
async def invoke(clisplit):
    from fractal.client.client import handle

    async def __invoke(args: str):
        return await handle(clisplit(args))

    return __invoke


@pytest.fixture
def clear_task_cache():
    from fractal.client.config import settings

    cache_dir = str(Path(f"{settings.FRACTAL_CACHE_PATH}").expanduser())
    cache_file = f"{cache_dir}/tasks"
    if os.path.isfile(cache_file):
        os.remove(cache_file)


from .fixtures_testserver import *  # noqa: 401
