import asyncio
import logging
from os import environ

import pytest


@pytest.fixture(scope="session")
def temp_data_dir(tmp_path_factory):
    yield tmp_path_factory.mktemp("data")


@pytest.fixture(scope="session")
def temp_db_path(tmp_path_factory):
    db_dir = tmp_path_factory.mktemp("db")
    yield db_dir / "test.db"


@pytest.fixture(scope="function")
def clear_db():
    """
    Clear the db without dropping the schema

    ref: https://stackoverflow.com/a/5003705/283972
    """
    import contextlib
    from fractal_server.app.db import DB
    from sqlmodel import SQLModel

    meta = SQLModel.metadata
    with contextlib.closing(DB.engine_sync().connect()) as conn:
        trans = conn.begin()
        for table in reversed(meta.sorted_tables):
            conn.execute(table.delete())
        trans.commit()


@pytest.fixture(scope="session")
async def testserver(temp_data_dir, temp_db_path):
    # cf. https://stackoverflow.com/a/57816608/283972
    environ["JWT_SECRET_KEY"] = "secret_key"
    environ["DEPLOYMENT_TYPE"] = "development"
    environ["DATA_DIR_ROOT"] = temp_data_dir.as_posix()
    environ["FRACTAL_ROOT"] = temp_data_dir.as_posix()
    environ["FRACTAL_LOGGING_LEVEL"] = str(logging.DEBUG)

    environ["DB_ENGINE"] = "sqlite"

    tmp_db_path = temp_db_path
    environ["SQLITE_PATH"] = tmp_db_path.as_posix()

    import uvicorn
    from fractal_server.main import start_application
    from multiprocessing import Process

    # INIT DB
    from fractal_server.app.db import DB
    from sqlmodel import SQLModel

    # import fractal_server.app.models  # noqa: F401

    SQLModel.metadata.create_all(DB.engine_sync())

    # We are explicitly calling start_application() to bypass the task
    # collection routine
    app = start_application()

    config = uvicorn.Config(app, port=10080, log_level="debug")
    server = uvicorn.Server(config)

    def run_server():
        asyncio.run(server.serve())

    proc = Process(target=run_server, args=(), daemon=True)
    proc.start()
    yield environ["FRACTAL_SERVER"]
    proc.kill()


@pytest.fixture()
async def user_factory(client, testserver):
    async def __register_user(email: str, password: str, slurm_user: str):
        res = await client.post(
            f"{testserver}/auth/register",
            json=dict(email=email, password=password, slurm_user=slurm_user),
        )
        from devtools import debug

        debug(res.json())
        assert res.status_code == 201
        return res.json()

    return __register_user


@pytest.fixture
async def register_user(user_factory):
    return await user_factory(
        email=environ["FRACTAL_USER"],
        password=environ["FRACTAL_PASSWORD"],
        slurm_user=environ["SLURM_USER"],
    )
