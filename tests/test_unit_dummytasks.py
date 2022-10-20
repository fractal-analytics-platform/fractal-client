import asyncio
import json

from devtools import debug

from fractal_server.tasks import dummy as dummy_module
from fractal_server.tasks.dummy import dummy


FIRST_TEST_MESSAGE = "first call"
SECOND_TEST_MESSAGE = "second call"


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
    assert metadata_update == {"dummy": "dummy 0"}
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
    assert metadata_update == {"dummy": "dummy 1"}
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
    assert metadata_update == {"dummy": "dummy 0"}
    with open(out_path, "r") as f:
        data = json.load(f)
    debug(data)
