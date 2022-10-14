import importlib
import os


def test_import_parsl_backend():
    import fractal_server.app.runner.parsl

    os.environ.pop("DEPLOYMENT_TYPE")
    importlib.reload(fractal_server.config)
    import fractal_server.app.runner.parsl
