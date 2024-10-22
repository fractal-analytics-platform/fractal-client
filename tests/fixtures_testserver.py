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


@pytest.fixture(scope="session")
def override_server_settings(tmp_path_factory):

    from fractal_server.config import Settings, get_settings
    from fractal_server.syringe import Inject

    settings = Settings()

    settings.DB_ENGINE = "postgres-psycopg"
    settings.POSTGRES_DB = "fractal_client_test"
    settings.POSTGRES_USER = "postgres"
    settings.POSTGRES_PASSWORD = "postgres"

    settings.FRACTAL_RUNNER_BACKEND = "local"

    settings.JWT_SECRET_KEY = "secret_key"
    base_folder = tmp_path_factory.mktemp("tmp_path")
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


@pytest.fixture(scope="session", autouse=True)
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
def task_factory(invoke):
    def _task_factory(
        name: str,
        command_non_parallel: Optional[str] = None,
        command_parallel: Optional[str] = None,
        version: Optional[str] = None,
        meta_non_parallel: Optional[str] = None,
        meta_parallel: Optional[str] = None,
        args_schema_non_parallel: Optional[str] = None,
        args_schema_parallel: Optional[str] = None,
        args_schema_version: Optional[str] = None,
    ):
        cmd = "task new"
        if command_non_parallel is not None:
            cmd += f" --command-non-parallel {command_non_parallel}"
        if command_parallel is not None:
            cmd += f" --command-parallel {command_parallel}"
        if version is not None:
            cmd += f" --version {version}"
        if meta_non_parallel is not None:
            cmd += f" --meta-non-parallel {meta_non_parallel}"
        if meta_parallel is not None:
            cmd += f" --meta-parallel {meta_parallel}"
        if args_schema_non_parallel is not None:
            cmd += f" --args-schema-non-parallel {args_schema_non_parallel}"
        if args_schema_parallel is not None:
            cmd += f" --args-schema-parallel {args_schema_parallel}"
        if args_schema_version is not None:
            cmd += f" --args-schema-version {args_schema_version}"
        cmd += f" {name}"

        res = invoke(cmd)
        return res.data

    return _task_factory


@pytest.fixture
def project_factory(invoke):
    def _project_factory(name: str):
        res = invoke(f"project new {name}")
        return res.data

    return _project_factory


@pytest.fixture
def workflow_factory(invoke):
    def _workflow_factory(name: str, project_id: int):
        res = invoke(f"worklow new {name} {project_id}")
        return res.data

    return _workflow_factory


@pytest.fixture
def job_factory(invoke):
    def _job_factory(
        project_id: int,
        workflow_id: int,
        dataset_id: int,
        start: Optional[int] = None,
        end: Optional[int] = None,
        worker_init: Optional[str] = None,
    ):
        cmd = "job submit"
        if start is not None:
            cmd += f" --start {start}"
        if end is not None:
            cmd += f" --end {end}"
        if worker_init is not None:
            cmd += f" --worker-init {worker_init}"
        cmd += f" {project_id} {workflow_id} {dataset_id}"

        res = invoke(cmd)
        return res.data

    return _job_factory


@pytest.fixture
def user_factory(client_superuser):
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
def register_user(user_factory, tmp_path, override_settings):
    user = tmp_path.as_posix().split("/")[-1]
    created_user = user_factory(
        email=f"{user}@fractal.xy", password=user, username=user
    )
    override_settings(FRACTAL_USER=f"{user}@fractal.xy", FRACTAL_PASSWORD=user)
    yield created_user
