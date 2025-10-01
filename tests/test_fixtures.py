from pathlib import Path

from devtools import debug
from fractal_client.config import settings


def test_clear_cache():
    cache_dir = settings.FRACTAL_CACHE_PATH
    debug(cache_dir)
    assert cache_dir != str(Path.home() / ".cache/fractal")


def test_override_settings(override_settings):
    override_settings(FRACTAL_CACHE_PATH="/tmp/xy")
    cache_dir = settings.FRACTAL_CACHE_PATH
    assert cache_dir == "/tmp/xy"
