from os import environ

import pytest
from asyncclick.testing import CliRunner
from httpx import AsyncClient


environ["FRACTAL_USER"] = "test@fake-exact-lab.it"
environ["FRACTAL_PASSWORD"] = "password"
environ["FRACTAL_SERVER"] = "http://127.0.0.1:10080"
environ["DB_ECHO"] = "0"


@pytest.fixture
async def cli():
    yield CliRunner()


@pytest.fixture
async def client():
    async with AsyncClient() as client:
        yield client


from .fixtures_testserver import *  # noqa: 401
