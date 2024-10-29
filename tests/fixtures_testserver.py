import logging
import os
import shlex
import subprocess
import time
from pathlib import Path
from typing import Optional

import pytest
from httpx import ConnectError

from fractal_client.client import handle

DB_NAME = "pytest-fractal-client"

logger = logging.getLogger("fractal-client")
logger.setLevel(logging.DEBUG)

PORT = 8765


@pytest.fixture
def superuser(invoke_as_superuser):
    return invoke_as_superuser("user whoami").data


@pytest.fixture(scope="session")
def tester():
    return dict(email="client_tester@example.org", password="pytest")


def _run_command(cmd: str) -> str:
    logging.warning(f"Now running {cmd=}")
    res = subprocess.run(
        shlex.split(cmd),
        capture_output=True,
        env=dict(PGPASSWORD="postgres", **os.environ),
        encoding="utf-8",
    )
    if res.returncode != 0:
        logging.error(f"{res.stdout=}")
        logging.error(f"{res.stderr=}")
        raise RuntimeError(res.stderr)
    else:
        return res.stdout


@pytest.fixture(scope="session", autouse=True)
def testserver(tester, tmpdir_factory, request):

    FRACTAL_TASK_DIR = str(tmpdir_factory.mktemp("TASKS"))
    FRACTAL_RUNNER_WORKING_BASE_DIR = str(tmpdir_factory.mktemp("JOBS"))

    env_file = Path(".fractal_server.env")
    with env_file.open("w") as f:
        f.write(
            "DB_ENGINE=postgres-psycopg\n"
            "POSTGRES_HOST=localhost\n"
            f"POSTGRES_DB={DB_NAME}\n"
            "POSTGRES_USER=postgres\n"
            "POSTGRES_PASSWORD=postgres\n"
            "FRACTAL_RUNNER_BACKEND=local\n"
            "JWT_SECRET_KEY=secret_key\n"
            f"FRACTAL_TASKS_DIR={FRACTAL_TASK_DIR}\n"
            "FRACTAL_RUNNER_WORKING_BASE_DIR="
            f"{FRACTAL_RUNNER_WORKING_BASE_DIR}\n"
            "FRACTAL_LOGGING_LEVEL=0\n"
        )
    _run_command(
        f"dropdb --username=postgres --host localhost --if-exists {DB_NAME}"
    )
    _run_command(f"createdb --username=postgres --host localhost {DB_NAME}")
    _run_command("poetry run fractalctl set-db")

    LOG_DIR = Path(
        os.environ.get(
            "GHA_FRACTAL_SERVER_LOG",
            tmpdir_factory.mktemp("LOGS"),
        ),
    )
    path_out = LOG_DIR / "server_out"
    path_err = LOG_DIR / "server_err"
    f_out = path_out.open("w")
    f_err = path_err.open("w")

    server_process = subprocess.Popen(
        shlex.split(f"poetry run fractalctl start --port {PORT}"),
        stdout=f_out,
        stderr=f_err,
    )

    # Wait until the server is up
    TIMEOUT = 8.0
    t_start = time.perf_counter()
    while True:
        try:
            res = handle(shlex.split("fractal version"))
            if res.retcode == 0:
                break
            else:
                raise ConnectError("fractal-server not ready")
        except ConnectError:
            logger.debug("Fractal server not ready, wait one more second.")
            if time.perf_counter() - t_start > TIMEOUT:
                raise RuntimeError(
                    f"Could not start up server within {TIMEOUT} seconds,"
                    " in `testserver` fixture."
                )
            time.sleep(0.1)

    handle(
        shlex.split(
            (
                "fractal --user admin@fractal.xy --password 1234 "
                f"user register {tester['email']} {tester['password']}"
            )
        )
    )

    yield

    request.session.warn(
        Warning(
            f"\n\nTerminating Fractal Server (PID: {server_process.pid}).\n"
            f"stdout -> {path_out}\n"
            f"stderr -> {path_err}\n"
        )
    )

    server_process.terminate()
    server_process.kill()
    _run_command(f"dropdb --username=postgres --host localhost {DB_NAME}")
    env_file.unlink()
    f_out.close()
    f_err.close()


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
        res = invoke(f"workflow new {name} {project_id}")
        return res.data

    return _workflow_factory


@pytest.fixture
def dataset_factory(invoke):
    def _dataset_factory(
        project_id: int,
        name: str,
        zarr_dir: str,
        filters: Optional[str] = None,
    ):
        cmd = "project add-dataset"
        if filters is not None:
            cmd += f" --filters {filters}"
        cmd += f" {project_id} {name} {zarr_dir}"

        res = invoke(cmd)
        return res.data

    return _dataset_factory


@pytest.fixture
def user_factory(invoke_as_superuser):
    def __user_factory(
        email: str,
        password: str,
        cache_dir: Optional[str] = None,
        slurm_user: Optional[str] = None,
        username: Optional[str] = None,
        superuser: bool = False,
    ):
        cmd = "user register"
        if cache_dir is not None:
            cmd += f" --cache-dir {cache_dir}"
        if slurm_user is not None:
            cmd += f" --slurm-user {slurm_user}"
        if username is not None:
            cmd += f" --username {username}"
        if superuser is True:
            cmd += " --superuser"
        cmd += f" {email} {password}"

        res = invoke_as_superuser(cmd)
        return res.data

    return __user_factory
