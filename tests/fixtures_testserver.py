import logging
import os
import signal
import subprocess
import time
from typing import Optional

import httpx
import pytest


logger = logging.getLogger("fractal-client")
logger.setLevel(logging.DEBUG)

PORT = 10080


@pytest.fixture(scope="session", autouse=True)
def testserver():

    server_process = subprocess.Popen(
        ["bash", "tests/run_test_server.sh"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

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

    try:
        yield server_process
    finally:
        if server_process.poll() is None:
            os.kill(server_process.pid, signal.SIGTERM)
            server_process.wait()


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
