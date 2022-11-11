import asyncio
import json
from pathlib import Path
from typing import Any
from typing import Dict
from typing import Optional

import pytest
from devtools import debug
from pydantic import BaseModel


class MockTask(BaseModel):
    name: str
    command: str
    parallelization_level: Optional[str] = None


class MockWorkflowTask(BaseModel):
    order: int = 0
    task: MockTask
    arguments: Dict = {}

    @property
    def is_parallel(self) -> bool:
        return bool(self.task.parallelization_level)

    @property
    def parallelization_level(self) -> Optional[str]:
        return self.task.parallelization_level

    def assemble_args(self, extra: Dict[str, Any] = None):
        """
        Merge of `extra` arguments and `self.arguments`.

        Return
        ------
        full_arsgs (Dict):
            A dictionary consisting of the merge of `extra` and
            self.arguments.
        """
        full_args = {}
        if extra:
            full_args.update(extra)
        full_args.update(self.arguments)
        return full_args


async def execute_command(cmd, **kwargs):
    proc = await asyncio.create_subprocess_shell(
        cmd,
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.PIPE,
        **kwargs,
    )
    stdout, stderr = await proc.communicate()
    debug(cmd, stdout, stderr)
    if proc.returncode != 0:
        raise RuntimeError(stderr.decode("UTF-8"))
    return stdout.decode("UTF-8").strip()


@pytest.fixture
async def dummy_task_package(testdata_path, tmp_path) -> Path:
    """
    Yields
    ------
    wheel_path : Path
        the path to the built wheel package
    """
    from fractal_server import tasks as task_package

    PACKAGE_TEMPLATE_PATH = testdata_path / "fractal-tasks-dummy"
    PACKAGE_PATH = tmp_path / "fractal-tasks-dummy"
    SOURCE_PATH = PACKAGE_PATH / "fractal_tasks_dummy"
    DUMMY_PACKAGE = Path(task_package.__file__).parent

    # copy template to temp
    await execute_command(f"cp -r {PACKAGE_TEMPLATE_PATH} {PACKAGE_PATH}")
    # copy content of task_package to PACKAGE_PATH
    await execute_command(f"cp {DUMMY_PACKAGE}/*.* {SOURCE_PATH}")
    await execute_command("poetry build", cwd=PACKAGE_PATH)
    wheel_relative = await execute_command("ls dist/*.whl", cwd=PACKAGE_PATH)
    wheel_path = PACKAGE_PATH / wheel_relative
    yield wheel_path


@pytest.fixture
async def dummy_task_package_invalid_manifest(testdata_path, tmp_path) -> Path:
    """
    Yields
    ------
    wheel_path : Path
        the path to the built wheel package
    """
    from fractal_server import tasks as task_package

    PACKAGE_TEMPLATE_PATH = testdata_path / "fractal-tasks-dummy"
    PACKAGE_PATH = tmp_path / "fractal-tasks-dummy"
    SOURCE_PATH = PACKAGE_PATH / "fractal_tasks_dummy"
    DUMMY_PACKAGE = Path(task_package.__file__).parent

    # copy template to temp
    await execute_command(f"cp -r {PACKAGE_TEMPLATE_PATH} {PACKAGE_PATH}")
    # copy content of task_package to PACKAGE_PATH
    await execute_command(f"cp {DUMMY_PACKAGE}/*.* {SOURCE_PATH}")

    # Make manifest invalid
    MANIFEST_PATH = f"{SOURCE_PATH}/__FRACTAL_MANIFEST__.json"
    with open(MANIFEST_PATH, "r") as fin:
        manifest = json.load(fin)
    task_list = manifest["task_list"]
    task_list[0].pop("default_args")
    manifest["task_list"] = task_list
    with open(MANIFEST_PATH, "w") as fout:
        json.dump(manifest, fout)

    await execute_command("poetry build", cwd=PACKAGE_PATH)
    wheel_relative = await execute_command("ls dist/*.whl", cwd=PACKAGE_PATH)
    wheel_path = PACKAGE_PATH / wheel_relative
    yield wheel_path
