import json
import logging
import os
import shlex
import subprocess
import sys
import time
from pathlib import Path

import pytest
from fractal_client.client import handle
from fractal_client.interface import Interface
from httpx import ConnectError

DB_NAME = "pytest-fractal-client"
FRACTAL_SERVER_PORT = 8765

logger = logging.getLogger("fractal-client")
logger.setLevel(logging.DEBUG)


def _run_command(cmd: str) -> str:
    logging.warning(f"Now running {cmd=}")
    env = os.environ
    if "PGPASSWORD" not in os.environ:
        env["PGPASSWORD"] = "postgres"
    res = subprocess.run(
        shlex.split(cmd),
        capture_output=True,
        env=env,
        encoding="utf-8",
    )
    if res.returncode != 0:
        logging.error(f"{res.stdout=}")
        logging.error(f"{res.stderr=}")
        raise RuntimeError(res.stderr)
    else:
        return res.stdout


def _drop_db():
    """
    Note: `--force` helps in case the `fractal-server` process did not
    terminate properly, thus leaving some open database connections.
    """
    _run_command(
        (
            "dropdb --username=postgres --host localhost "
            f"--if-exists {DB_NAME} --force"
        )
    )


def _split_and_handle(cli_string: str) -> Interface:
    return handle(shlex.split(cli_string))


@pytest.fixture(scope="session", autouse=True)
def testserver(tester, tmpdir_factory, request):
    FRACTAL_TASK_DIR = str(tmpdir_factory.mktemp("TASKS"))
    FRACTAL_RUNNER_WORKING_BASE_DIR = str(tmpdir_factory.mktemp("JOBS"))
    ADMIN_EMAIL = "admin@fractal.xy"
    ADMIN_PWD = "1234"
    ADMIN_PROJECT_DIR = str(tmpdir_factory.mktemp("ADMIN_PROJECT"))

    env_file = Path(".fractal_server.env")
    with env_file.open("w") as f:
        f.write(
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
            "FRACTAL_VIEWER_AUTHORIZATION_SCHEME=viewer-paths\n"
        )
    _drop_db()
    _run_command(f"createdb --username=postgres --host localhost {DB_NAME}")
    _run_command("uv run fractalctl set-db")
    _run_command(
        "uv run fractalctl init-db-data "
        "--resource default "
        "--profile default "
        f"--admin-email {ADMIN_EMAIL} "
        f"--admin-pwd {ADMIN_PWD} "
        f"--admin-project-dir {ADMIN_PROJECT_DIR}"
    )

    LOG_DIR = os.environ.get(
        "GHA_FRACTAL_SERVER_LOG",
        tmpdir_factory.mktemp("LOGS"),
    )

    path_out = Path(LOG_DIR, "server_out")
    path_err = Path(LOG_DIR, "server_err")
    f_out = path_out.open("w")
    f_err = path_err.open("w")

    server_process = subprocess.Popen(
        shlex.split(f"uv run fractalctl start --port {FRACTAL_SERVER_PORT}"),
        stdout=f_out,
        stderr=f_err,
    )

    # Wait until the server is up
    TIMEOUT = 8.0
    t_start = time.perf_counter()
    while True:
        try:
            res = _split_and_handle("fractal version")
            if "refused" not in res.data:
                break
            else:
                f_out.close()
                f_err.close()
                raise ConnectError("fractal-server not ready")
        except ConnectError:
            logger.debug("Fractal server not ready, wait one more second.")
            if time.perf_counter() - t_start > TIMEOUT:
                raise RuntimeError(
                    f"Could not start up server within {TIMEOUT} seconds,"
                    " in `testserver` fixture."
                )
            time.sleep(0.1)

    res = _split_and_handle(
        "fractal --user admin@fractal.xy --password 1234 --batch "
        "user register "
        f"{tester['email']} {tester['password']} {tester['project_dir']}"
    )
    user_id = res.data

    current_py_version = f"{sys.version_info.major}.{sys.version_info.minor}"
    RES_PROF_DIR = str(tmpdir_factory.mktemp("resource_and_profile"))

    resource_json = Path(RES_PROF_DIR) / "resource.json"
    profile_json = Path(RES_PROF_DIR) / "profile.json"

    resource = dict(
        name="local resource tester",
        type="local",
        jobs_local_dir=(Path(RES_PROF_DIR) / "jobs").as_posix(),
        tasks_local_dir=(Path(RES_PROF_DIR) / "tasks").as_posix(),
        jobs_runner_config={"parallel_tasks_per_job": 1},
        tasks_python_config={
            "default_version": current_py_version,
            "versions": {
                current_py_version: sys.executable,
            },
        },
        tasks_pixi_config={},
        jobs_poll_interval=0,
    )
    with Path(resource_json).open("w") as f:
        json.dump(resource, f)

    res = _split_and_handle(
        "fractal --user admin@fractal.xy --password 1234 --batch resource new "
        f"{resource_json}"
    )
    resource_id = res.data

    profile = dict(
        name="local_resource_profile_objects",
        resource_id=resource_id,
        resource_type="local",
    )
    with Path(profile_json).open("w") as f:
        json.dump(profile, f)

    res = _split_and_handle(
        "fractal --user admin@fractal.xy --password 1234 --batch profile new "
        f"{resource_id} {profile_json}"
    )
    profile_id = res.data

    _split_and_handle(
        "fractal --user admin@fractal.xy --password 1234 user edit "
        f"--new-profile-id {profile_id} {user_id}"
    )

    yield

    request.session.warn(
        Warning(
            f"\n\nTerminating Fractal Server (PID: {server_process.pid}).\n"
            f"stdout -> {path_out}\n"
            f"stderr -> {path_err}\n"
        )
    )

    # Start cleanup

    server_process.poll()
    server_process.terminate()
    server_process.kill()
    server_process.poll()

    f_out.close()
    f_err.close()
    env_file.unlink()

    _drop_db()


@pytest.fixture
def task_factory(invoke):
    def _task_factory(
        name: str,
        command_non_parallel: str | None = None,
        command_parallel: str | None = None,
        version: str | None = None,
        meta_non_parallel: str | None = None,
        meta_parallel: str | None = None,
        args_schema_non_parallel: str | None = None,
        args_schema_parallel: str | None = None,
        args_schema_version: str | None = None,
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
        if command_non_parallel is not None and command_parallel is not None:
            cmd += " --task-type compound"
        elif command_non_parallel is not None:
            cmd += " --task-type non_parallel"
        else:
            cmd += " --task-type parallel"

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
    ):
        cmd = "project add-dataset"
        cmd += f" {project_id} {name} --zarr-dir {zarr_dir}"

        res = invoke(cmd)
        return res.data

    return _dataset_factory


@pytest.fixture
def user_factory(invoke_as_superuser):
    def __user_factory(
        email: str,
        password: str,
        project_dir: str = "/fake",
        username: str | None = None,
        superuser: bool = False,
    ):
        cmd = "user register"
        if username is not None:
            cmd += f" --username {username}"
        if superuser is True:
            cmd += " --superuser"
        cmd += f" {email} {password} {project_dir}"

        res = invoke_as_superuser(cmd)
        return res.data

    return __user_factory
