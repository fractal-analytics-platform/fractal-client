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


@pytest.fixture
async def db(testserver):
    """
    NOTE: Only use this fixture within other fixtures!!!
    """
    from fractal_server.app.db import get_db

    async for db in get_db():
        yield db


@pytest.fixture
async def task_factory(db):
    from fractal_server.app.models.task import Task

    async def _task_factory(**task_args_override):
        task_args = dict(
            name="test_task",
            command="cmd",
            source="source",
            input_type="Any",
            output_type="Any",
        )
        task_args.update(task_args_override)
        t = Task(**task_args)
        db.add(t)
        await db.commit()
        await db.refresh(t)
        return t

    return _task_factory


@pytest.fixture
async def project_factory(db):
    from fractal_server.app.models.project import Project

    async def _project_factory(user_id=None, **project_args_override):
        project_args = dict(name="name", project_dir="project/dir")
        project_args.update(project_args_override)
        p = Project(**project_args)
        if user_id:
            from fractal_server.app.security import User

            user = await db.get(User, user_id)
            p.user_member_list.append(user)
        db.add(p)
        await db.commit()
        await db.refresh(p)
        return p

    return _project_factory


@pytest.fixture
async def workflow_factory(db, project_factory, register_user):
    from fractal_server.app.models.workflow import Workflow

    async def _workflow_factory(**wf_args_override):
        if "project_id" not in wf_args_override:
            p = await project_factory(user_id=register_user["id"])
            wf_args_override["project_id"] = p.id

        wf_args = dict(
            name="name",
        )
        wf_args.update(wf_args_override)
        wf = Workflow(**wf_args)
        db.add(wf)
        await db.commit()
        await db.refresh(wf)
        return wf

    return _workflow_factory


@pytest.fixture()
async def user_factory(client, testserver):
    async def __register_user(email: str, password: str, slurm_user: str):
        res = await client.post(
            f"{testserver}/auth/register",
            json=dict(email=email, password=password, slurm_user=slurm_user),
        )
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
