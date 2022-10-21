import os
import sys

import pytest

try:
    import parsl  # noqa: F401

    HAS_PARSL = True
except ImportError:
    HAS_PARSL = False


@pytest.mark.skipif(
    not HAS_PARSL, reason="Optional dependency `Parsl` is not installed"
)
def test_import_parsl_backend():
    import fractal_server.config

    os.environ.pop("DEPLOYMENT_TYPE")
    del sys.modules["fractal_server.config"]
    del fractal_server.config

    import fractal_server.app.runner.parsl
