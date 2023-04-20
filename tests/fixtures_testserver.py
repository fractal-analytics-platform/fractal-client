import asyncio
import logging
from os import environ
from typing import Optional

import pytest


logger = logging.getLogger("fractal-client")
logger.setLevel(logging.DEBUG)


@pytest.fixture
def override_server_settings(tmp_path):
    from fractal_server.config import Settings, get_settings
    from fractal_server.syringe import Inject

    settings = Settings()

    tmp_db_path = tmp_path / "db/test.db"
    tmp_db_path.parent.mkdir()
    settings.DB_ENGINE = "sqlite"
    settings.SQLITE_PATH = tmp_db_path

    settings.JWT_SECRET_KEY = "secret_key"
    settings.DEPLOYMENT_TYPE = "development"
    base_folder = tmp_path
    settings.FRACTAL_TASKS_DIR = base_folder / "FRACTAL_TASKS_DIR"
    settings.FRACTAL_RUNNER_WORKING_BASE_DIR = (
        base_folder / "FRACTAL_RUNNER_WORKING_BASE_DIR"
    )
    settings.FRACTAL_LOGGING_LEVEL = logging.DEBUG

    def _get_settings():
        return settings

    Inject.override(get_settings, _get_settings)
    try:
        yield
    finally:
        Inject.pop(get_settings)


@pytest.fixture(scope="function", autouse=True)
async def testserver(override_server_settings):
    import uvicorn
    from fractal_server.main import start_application
    from fractal_server.main import _create_first_user
    from multiprocessing import Process
    from fractal_server.app.db import DB
    from fractal_server.app.models import SQLModel

    # INIT DB
    DB.set_db()
    logger.debug(DB.engine_sync().url)
    SQLModel.metadata.create_all(DB.engine_sync())

    # We are explicitly calling start_application() to bypass the task
    # collection routine
    app = start_application()
    await _create_first_user("admin@fractal.xy", "1234", is_superuser=True)

    config = uvicorn.Config(app, port=10080, log_level="debug")
    server = uvicorn.Server(config)

    def run_server():
        asyncio.run(server.serve())

    # Running testserver in a separate process
    # cf. https://stackoverflow.com/a/57816608/283972
    proc = Process(target=run_server, args=(), daemon=True)
    proc.start()
    import time

    time.sleep(0.2)  # NOTE required to let the server comes up
    logger.debug(environ["FRACTAL_SERVER"])
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
        project_args = dict(name="name")
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


@pytest.fixture
async def job_factory(db):
    from fractal_server.app.models.job import ApplyWorkflow

    async def _job_factory(**job_args_override):
        job_args = dict(
            project_id=1,
            input_dataset_id=1,
            output_dataset_id=2,
            workflow_id=1,
            overwrite_input=False,
            worker_init="WORKER_INIT string",
        )
        job_args.update(job_args_override)
        j = ApplyWorkflow(**job_args)
        db.add(j)
        await db.commit()
        await db.refresh(j)
        return j

    return _job_factory


@pytest.fixture
async def user_factory(client_superuser, testserver):
    async def __register_user(
        email: str, password: str, slurm_user: Optional[str] = None
    ):
        payload = dict(email=email, password=password)
        if slurm_user:
            payload["slurm_user"] = slurm_user
        res = await client_superuser.post(
            f"{testserver}/auth/register",
            json=payload,
        )
        assert res.status_code == 201
        return res.json()

    return __register_user


@pytest.fixture
async def register_user(user_factory):
    return await user_factory(
        email=environ["FRACTAL_USER"],
        password=environ["FRACTAL_PASSWORD"],
    )
