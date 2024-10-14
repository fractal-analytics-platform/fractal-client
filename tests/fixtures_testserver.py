import logging
import time
from multiprocessing import Process
from os import environ
from typing import Optional

import httpx
import pytest
import uvicorn
from sqlmodel import select


logger = logging.getLogger("fractal-client")
logger.setLevel(logging.DEBUG)

PORT = 10080


@pytest.fixture
def override_server_settings(tmp_path):
    from fractal_server.config import Settings, get_settings
    from fractal_server.syringe import Inject

    settings = Settings()

    settings.DB_ENGINE = "postgres-psycopg"
    settings.POSTGRES_DB = "fractal_client_test"
    settings.POSTGRES_USER = "postgres"
    settings.POSTGRES_PASSWORD = "postgres"

    settings.FRACTAL_RUNNER_BACKEND = "local"

    settings.JWT_SECRET_KEY = "secret_key"
    base_folder = tmp_path
    settings.FRACTAL_TASKS_DIR = base_folder / "FRACTAL_TASKS_DIR"
    settings.FRACTAL_RUNNER_WORKING_BASE_DIR = (
        base_folder / "FRACTAL_RUNNER_WORKING_BASE_DIR"
    )
    settings.FRACTAL_LOGGING_LEVEL = logging.DEBUG
    settings.FRACTAL_API_SUBMIT_RATE_LIMIT = 0

    def _get_settings():
        return settings

    Inject.override(get_settings, _get_settings)
    try:
        yield
    finally:
        Inject.pop(get_settings)


@pytest.fixture(scope="function", autouse=True)
def testserver(override_server_settings):

    from fractal_server.app.db import DB
    from fractal_server.app.models.security import SQLModel
    from fractal_server.app.models.security import UserOAuth
    from fractal_server.app.models.user_settings import UserSettings
    from fractal_server.app.models.security import UserGroup
    from fractal_server.app.models.linkusergroup import LinkUserGroup
    from fractal_server.app.security import _create_first_group

    # INIT DB
    engine_sync = DB.engine_sync()
    logger.debug(engine_sync.url)
    SQLModel.metadata.create_all(engine_sync)

    # Create default group and first superuser
    # NOTE: we have to do it here, because we are not calling the `set_db`
    # function from fractal-server. This would change with
    # https://github.com/fractal-analytics-platform/fractal-client/issues/697
    # NOTE: `hashed_password` is the bcrypt hash of "1234", see
    # https://github.com/fractal-analytics-platform/fractal-server/issues/1750
    _create_first_group()
    with next(DB.get_sync_db()) as db:
        user = UserOAuth(
            email="admin@fractal.xy",
            hashed_password=(
                "$2b$12$K0C4t7XILgpcQx35V3QE3enOODQ1IH9pzW49nqjHbrx2uQTMVYsQC"
            ),
            username=environ["FRACTAL_USERNAME"],
            is_superuser=True,
            is_verified=True,
            is_active=True,
        )
        empty_user_settings = UserSettings()
        user.settings = empty_user_settings
        db.add(user)
        db.commit()

        first_group = db.execute(select(UserGroup)).scalar()
        first_user = db.execute(select(UserOAuth)).scalar()

        link = LinkUserGroup(group_id=first_group.id, user_id=first_user.id)
        db.add(link)
        db.commit()

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

    # Cleanup DB
    engine_sync.dispose()
    try:
        DB._engine_async
        raise
    except AttributeError:
        # we show here that we do not need to dispose of `engine_async`,
        # because it is never used.
        pass
    SQLModel.metadata.drop_all(engine_sync)
    logger.debug("Dropped all tables from the database.")

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
    from fractal_server.app.models.v2.task import TaskV2
    from fractal_server.app.models.v2.task import TaskGroupV2

    def _task_factory(user_id: int, **task_args_override):
        task_args = dict(name="test_task", type="parallel")
        task_args.update(task_args_override)
        t = TaskV2(**task_args)

        db.add(
            TaskGroupV2(
                user_id=user_id,
                origin="other",
                pkg_name=t.name,
                task_list=[t],
            )
        )

        db.commit()
        db.refresh(t)
        return t

    return _task_factory


@pytest.fixture
def project_factory(db):
    from fractal_server.app.models.v2.project import ProjectV2

    def _project_factory(user_id=None, **project_args_override):
        project_args = dict(name="name")
        project_args.update(project_args_override)
        p = ProjectV2(**project_args)
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
    from fractal_server.app.models.v2.workflow import WorkflowV2

    def _workflow_factory(**wf_args_override):
        wf_args = dict(name="name")
        wf_args.update(wf_args_override)
        wf = WorkflowV2(**wf_args)
        db.add(wf)
        db.commit()
        db.refresh(wf)
        return wf

    return _workflow_factory


@pytest.fixture
def job_factory(db):
    from fractal_server.app.models.v2.job import JobV2
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
            dataset_dump=dict(
                id=1,
                name="ds-in",
                zarr_dir="/abc",
                project_id=1,
                timestamp_created=str(get_timestamp()),
                filters=dict(attributes=dict(a=1), types=dict(b=True)),
            ),
            project_dump=dict(
                id=1,
                name="proj",
                timestamp_created=str(get_timestamp()),
            ),
            start_timestamp=get_timestamp(),
            user_email="test@test.test",
        )
        job_args.update(job_args_override)
        j = JobV2(**job_args)
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
def register_user(user_factory):

    created_user = user_factory(
        email=environ["FRACTAL_USER"],
        password=environ["FRACTAL_PASSWORD"],
        username=environ["FRACTAL_USERNAME"],
    )

    yield created_user
