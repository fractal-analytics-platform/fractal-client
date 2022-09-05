import pytest
from asyncclick.testing import CliRunner


@pytest.fixture
async def cli():
    yield CliRunner()
