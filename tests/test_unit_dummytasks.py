import asyncio
import json
from pathlib import Path

import pytest
from devtools import debug

from fractal_server import tasks as tasks_package
from fractal_server.app.schemas.manifest import ManifestV1
from fractal_server.tasks import dummy as dummy_module
from fractal_server.tasks.dummy import dummy
from fractal_server.tasks.dummy_parallel import dummy_parallel


FIRST_TEST_MESSAGE = "first call"
SECOND_TEST_MESSAGE = "second call"
ERROR_MESSAGE = "This is an error message"


def test_dummy_direct_call(tmp_path):
    out_path = tmp_path / "out.json"
    metadata_update = dummy(
        input_paths=[tmp_path],
        output_path=out_path,
        metadata={"before": "test"},
        message=FIRST_TEST_MESSAGE,
        index=0,
    )
    assert out_path.exists()
    assert metadata_update == {"dummy": "dummy 0", "index": [0, 1, 2]}
    with open(out_path, "r") as f:
        data = json.load(f)
    debug(data)

    assert len(data) == 1
    assert data[0]["message"] == FIRST_TEST_MESSAGE

    # Second call
    metadata_update = dummy(
        input_paths=[tmp_path],
        output_path=out_path,
        metadata={"before": "test"},
        message=SECOND_TEST_MESSAGE,
        index=1,
    )
    assert metadata_update == {"dummy": "dummy 1", "index": [0, 1, 2]}
    with open(out_path, "r") as f:
        data = json.load(f)
    debug(data)

    assert len(data) == 2
    assert data[1]["message"] == SECOND_TEST_MESSAGE


async def test_dummy_process_call(tmp_path):
    args_file = tmp_path / "args.json"
    out_path = tmp_path / "out.json"
    args = dict(
        input_paths=[tmp_path.as_posix()],
        output_path=out_path.as_posix(),
        metadata={"before": "test"},
        message=FIRST_TEST_MESSAGE,
        index=0,
    )
    with open(args_file, "w") as fargs:
        json.dump(args, fargs)

    cmd = f"python {dummy_module.__file__} -j {args_file}"
    proc = await asyncio.create_subprocess_shell(
        cmd,
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.PIPE,
    )
    stdout, stderr = await proc.communicate()
    debug(proc.returncode)
    debug(stdout)
    debug(stderr)
    metadata_update = json.loads(stdout)
    assert proc.returncode == 0

    assert out_path.exists()
    assert metadata_update == {"dummy": "dummy 0", "index": [0, 1, 2]}
    with open(out_path, "r") as f:
        data = json.load(f)
    debug(data)


def test_dummy_fail_direct_call(tmp_path):
    out_path = tmp_path / "out.json"
    with pytest.raises(ValueError):
        dummy(
            input_paths=[tmp_path],
            output_path=out_path,
            metadata={"before": "test"},
            message=ERROR_MESSAGE,
            raise_error=True,
        )


async def test_dummy_fail_process_call(tmp_path):
    args_file = tmp_path / "args.json"
    out_path = tmp_path / "out.json"
    args = dict(
        input_paths=[tmp_path.as_posix()],
        output_path=out_path.as_posix(),
        metadata={"before": "test"},
        message=ERROR_MESSAGE,
        raise_error=True,
    )
    with open(args_file, "w") as fargs:
        json.dump(args, fargs)

    cmd = f"python {dummy_module.__file__} -j {args_file}"
    proc = await asyncio.create_subprocess_shell(
        cmd,
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.PIPE,
    )
    stdout, stderr = await proc.communicate()
    debug(proc.returncode)
    debug(stdout)
    debug(stderr)
    assert proc.returncode == 1
    assert f"ValueError: {ERROR_MESSAGE}" in str(stderr)


def test_dummy_parallel_direct_call(tmp_path):
    list_components = ["A", "B", "C"]
    out_path = tmp_path / "*.json"

    for component in list_components:
        metadata_update = dummy_parallel(
            input_paths=[tmp_path],
            output_path=out_path,
            component=component,
            metadata={"before": "test"},
            message=FIRST_TEST_MESSAGE,
        )
        assert not metadata_update

    assert out_path.parent.exists()
    out_files = list(out_path.parent.glob("*"))
    debug(out_files)
    assert len(out_files) == len(list_components)

    for out_file in out_files:
        with out_file.open("r") as fin:
            data = json.load(fin)
        assert out_file.name == f'{data["component"]}.json'
        assert data["message"] == FIRST_TEST_MESSAGE


def test_dummy_parallel_fail_direct_call(tmp_path):
    list_components = ["A", "B", "C"]
    out_path = tmp_path / "*.json"

    for component in list_components:
        with pytest.raises(ValueError):
            dummy_parallel(
                input_paths=[tmp_path],
                output_path=out_path,
                component=component,
                metadata={"before": "test"},
                message=ERROR_MESSAGE,
                raise_error=True,
            )


def test_manifest_validation():
    manifest_path = (
        Path(tasks_package.__file__).parent / "__FRACTAL_MANIFEST__.json"
    )
    with manifest_path.open("r") as f:
        manifest_dict = json.load(f)

    if manifest_dict["manifest_version"] == 1:
        manifest_obj = ManifestV1(**manifest_dict)
    debug(manifest_obj)
