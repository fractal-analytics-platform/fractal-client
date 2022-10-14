import os
import sys


def test_import_parsl_backend():
    import fractal_server.config

    os.environ.pop("DEPLOYMENT_TYPE")
    del sys.modules["fractal_server.config"]
    del fractal_server.config

    import fractal_server.app.runner.parsl
