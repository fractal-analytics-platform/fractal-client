from pathlib import Path

import pytest


@pytest.fixture(scope="session")
async def testdata_path() -> Path:
    TEST_DIR = Path(__file__).parent
    return TEST_DIR / "data/"


from .fixtures_server import *  # noqa F403
from .fixtures_server import override_environment  # noqa E402


# override_environment(Path(__file__).parent / "data/")
