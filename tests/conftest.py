from pathlib import Path

import pytest


@pytest.fixture(scope="session")
async def testdata_path() -> Path:
    TEST_DIR = Path(__file__).parent
    return TEST_DIR / "data/"


from .fixtures_server import *  # noqa F403
from .fixtures_server import override_enironment  # noqa E402


override_enironment(Path(__file__).parent / "data/")
