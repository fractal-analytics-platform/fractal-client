import logging
from os import environ
from typing import Optional

import httpx
import pytest


logger = logging.getLogger("fractal-client")
logger.setLevel(logging.DEBUG)

PORT = 10080


@pytest.fixture
def override_server_settings(tmp_path):
    from fractal_server.config import Settings, get_settings
    from fractal_server.syringe import Inject

    settings = Settings()

    tmp_db_path = tmp_path / "db/test.db"
    tmp_db_path.parent.mkdir()
    settings.DB_ENGINE = "sqlite"
    settings.SQLITE_PATH = tmp_db_path
    settings.FRACTAL_RUNNER_BACKEND = "local"

    settings.JWT_SECRET_KEY = "secret_key"
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
def testserver(override_server_settings):
    import uvicorn
    from multiprocessing import Process
    from fractal_server.app.db import DB
    import time
    from fractal_server.app.models import SQLModel

    # INIT DB
    DB.set_sync_db()
    logger.debug(DB.engine_sync().url)
    SQLModel.metadata.create_all(DB.engine_sync())

    # Run testserver in a separate process
    # cf. https://stackoverflow.com/a/57816608/283972

    def run_server():
        uvicorn.run(
            "fractal_server.main:app",
            port=PORT,
            log_level="debug",
            timeout_keep_alive=10,
        )

    proc = Process(target=run_server, args=(), daemon=True)
    proc.start()

    # Wait until the server is up
    TIMEOUT = 8
    time_used = 0
    while True:
        try:
            res = httpx.get(f"http://localhost:{PORT}/api/alive/")
            assert res.status_code == 200
            break
        except httpx.ConnectError:
            logger.debug("Fractal server not ready, wait one more second.")
            time.sleep(1)
            time_used += 1
            if time_used > TIMEOUT:
                raise RuntimeError(
                    f"Could not start up server within {TIMEOUT} seconds,"
                    " in `testserver` fixture."
                )

    logger.debug(environ["FRACTAL_SERVER"])
    yield environ["FRACTAL_SERVER"]
    proc.kill()


@pytest.fixture
def db(testserver):
    """
    NOTE: Only use this fixture within other fixtures!!!
    """
    from fractal_server.app.db import get_sync_db

    for db in get_sync_db():
        yield db


@pytest.fixture
def task_factory(db):
    from fractal_server.app.models.task import Task

    def _task_factory(**task_args_override):
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
        db.commit()
        db.refresh(t)
        return t

    return _task_factory


@pytest.fixture
def project_factory(db):
    from fractal_server.app.models.project import Project

    def _project_factory(user_id=None, **project_args_override):
        project_args = dict(name="name")
        project_args.update(project_args_override)
        p = Project(**project_args)
        if user_id:
            from fractal_server.app.security import User

            user = db.get(User, user_id)
            p.user_list.append(user)
        db.add(p)
        db.commit()
        db.refresh(p)
        return p

    return _project_factory


@pytest.fixture
def workflow_factory(db, project_factory):
    from fractal_server.app.models.workflow import Workflow

    def _workflow_factory(**wf_args_override):
        wf_args = dict(name="name")
        wf_args.update(wf_args_override)
        wf = Workflow(**wf_args)
        db.add(wf)
        db.commit()
        db.refresh(wf)
        return wf

    return _workflow_factory


@pytest.fixture
def job_factory(db):
    from fractal_server.app.models.job import ApplyWorkflow
    from fractal_server.utils import get_timestamp

    def _job_factory(**job_args_override):
        job_args = dict(
            project_id=1,
            input_dataset_id=1,
            output_dataset_id=2,
            workflow_id=1,
            worker_init="WORKER_INIT string",
            first_task_index=9999,
            last_task_index=9999,
            workflow_dump={},
            input_dataset_dump=dict(
                id=1,
                name="ds-in",
                read_only=False,
                project_id=1,
                resource_list=[dict(path="/tmp", id=1, dataset_id=1)],
                timestamp_created=str(get_timestamp()),
            ),
            output_dataset_dump=dict(
                id=2,
                name="ds-out",
                read_only=False,
                project_id=1,
                resource_list=[dict(path="/tmp", id=1, dataset_id=2)],
                timestamp_created=str(get_timestamp()),
            ),
            project_dump=dict(
                id=1,
                name="proj",
                read_only=True,
                timestamp_created=str(get_timestamp()),
            ),
            start_timestamp=get_timestamp(),
            user_email="test@test.test",
        )
        job_args.update(job_args_override)
        j = ApplyWorkflow(**job_args)
        db.add(j)
        db.commit()
        db.refresh(j)
        return j

    return _job_factory


@pytest.fixture
def user_factory(testserver, db, client_superuser):
    def __register_user(
        email: str,
        password: str,
        slurm_user: Optional[str] = None,
        username: Optional[str] = None,
    ):
        # Prepare payload
        new_user = dict(email=email, password=password)
        if slurm_user:
            new_user["slurm_user"] = slurm_user
        if username:
            new_user["username"] = username
        # Register user via API call
        res = client_superuser.post(
            f"http://localhost:{PORT}/auth/register/",
            json=new_user,
        )
        assert res.status_code == 201
        user_id = res.json()["id"]
        # Make user verified via API call
        res = client_superuser.patch(
            f"http://localhost:{PORT}/auth/users/{user_id}/",
            json=dict(is_verified=True),
        )
        assert res.status_code == 200
        return res.json()

    return __register_user


@pytest.fixture
def register_user(user_factory, db):
    from fractal_server.app.models import UserOAuth

    created_user = user_factory(
        email=environ["FRACTAL_USER"],
        password=environ["FRACTAL_PASSWORD"],
        username=environ["FRACTAL_USERNAME"],
    )

    yield created_user

    db_user = db.get(UserOAuth, created_user["id"])
    db.delete(db_user)
    db.commit()
