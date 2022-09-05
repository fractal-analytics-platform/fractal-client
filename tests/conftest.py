import asyncio
import pytest
from asyncclick.testing import CliRunner
from multiprocessing import Process
from os import environ
from devtools import debug
from httpx import AsyncClient


@pytest.fixture
async def cli():
    yield CliRunner()


@pytest.fixture
async def testserver(tmp_path):
    # cf. https://stackoverflow.com/a/57816608/283972
    import uvicorn
    from fractal_server import start_application
    from multiprocessing import Process

    environ["JWT_SECRET_KEY"] = "secret_key"
    environ["DEPLOYMENT_TYPE"] = "development"
    environ["DATA_DIR_ROOT"] = tmp_path.as_posix()

    environ["DB_ENGINE"] = "sqlite"
    # Shared in memory database,
    # c.f., https://stackoverflow.com/a/38089822/283972
    environ["SQLITE_PATH"] = "file:cachedb?mode=memory&cache=shared"

    # We are explicitly calling start_application() to bypass the task
    # collection routine
    app = start_application()

    config = uvicorn.Config(app, port=10080, log_level="debug")
    server = uvicorn.Server(config)

    def run_server():
        asyncio.run(server.serve())

    proc = Process(target=run_server, args=(), daemon=True)
    proc.start()
    yield "http://localhost:10080"
    proc.kill()


@pytest.fixture
async def client():
    async with AsyncClient() as client:
        yield client


@pytest.fixture
async def user_factory(client, testserver):
    res = await client.post(
        f"{testserver}/auth/register",
        json=dict(email="me@exact-lab.it", password="password")
    )
    debug(res)
    assert res.status_code == 201
