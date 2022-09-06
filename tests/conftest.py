from os import environ

import pytest
from asyncclick.testing import CliRunner
from httpx import AsyncClient


environ["FRACTAL_USER"] = "testuser"
environ["FRACTAL_PASSWORD"] = "password"


@pytest.fixture
async def cli():
    yield CliRunner()


@pytest.fixture
async def client():
    async with AsyncClient() as client:
        yield client


from .fixtures_testserver import *  # noqa: 401
