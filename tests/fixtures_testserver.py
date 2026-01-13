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


def _resource_and_profile_ids(base_path: Path, resource_name: str):
    base_path.mkdir(parents=True, exist_ok=True)
    current_py_version = f"{sys.version_info.major}.{sys.version_info.minor}"

    resource_json = base_path / "resource.json"
    profile_json = base_path / "profile.json"

    resource = dict(
        name=resource_name,
        type="local",
        jobs_local_dir=(base_path / "jobs").as_posix(),
        tasks_local_dir=(base_path / "tasks").as_posix(),
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

    with resource_json.open("w") as f:
        json.dump(resource, f)

    res = _split_and_handle(
        "fractal --user admin@fractal.xy --password 1234 --batch resource new "
        f"{resource_json.as_posix()}"
    )
    resource_id = res.data

    profile = dict(
        name=f"local_resource_{resource_id}_profile_objects",
        resource_id=resource_id,
        resource_type="local",
    )

    with profile_json.open("w") as f:
        json.dump(profile, f)

    res = _split_and_handle(
        "fractal --user admin@fractal.xy --password 1234 --batch profile new "
        f"{resource_id} {profile_json.as_posix()}"
    )
    profile_id = res.data

    return resource_id, profile_id


@pytest.fixture(scope="session", autouse=True)
def testserver(tester, tmpdir_factory, request):
    # Wait until the server is up
    TIMEOUT = 8.0
    t_start = time.perf_counter()
    while True:
        try:
            res = _split_and_handle("fractal version")
            if "refused" not in res.data:
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

    res = _split_and_handle(
        "fractal --user admin@fractal.xy --password 1234 --batch "
        "user register "
        f"{tester['email']} {tester['password']} {tester['project_dir']}"
    )
    user_id = res.data
    _, profile_id = _resource_and_profile_ids(
        base_path=Path(tmpdir_factory.mktemp("resource-and-profile")),
        resource_name=f"resource-{user_id}",
    )

    _split_and_handle(
        "fractal --user admin@fractal.xy --password 1234 user edit "
        f"--new-profile-id {profile_id} {user_id}"
    )

    yield

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
def user_factory(invoke_as_superuser, tmp_path):
    def __user_factory(
        email: str,
        password: str,
        project_dir: str = "/fake",
        superuser: bool = False,
    ):
        cmd = "user register"
        if superuser is True:
            cmd += " --superuser"
        cmd += f" {email} {password} {project_dir}"

        res = invoke_as_superuser(cmd)
        user_id = res.data["id"]

        _, profile_id = _resource_and_profile_ids(
            base_path=tmp_path / "resource-and-profile",
            resource_name=f"resource-{user_id}",
        )

        invoke_as_superuser(
            f"user edit --new-profile-id {profile_id} {user_id}"
        )

        res = invoke_as_superuser(f"user show {user_id}")
        return res.data

    return __user_factory
