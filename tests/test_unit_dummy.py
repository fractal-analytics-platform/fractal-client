import json
from pathlib import Path

from devtools import debug

from fractal.tasks.dummy import dummy


def test_dummy_task(tmp_path):
    input_path = [Path("/fake/input/path")]
    output_path = tmp_path
    outfile = output_path / "0.json"

    FIRST_MESSAGE = "first run"
    metadata = dummy(
        input_paths=input_path,
        output_path=output_path,
        message=FIRST_MESSAGE,
        metadata=None,
    )
    debug(metadata)

    with open(outfile, "r") as f:
        debug(f.read())
        f.seek(0)
        data = json.load(f)
    assert len(data) == 1
    assert data[0]["message"] == FIRST_MESSAGE

    SECOND_MESSAGE = "second run"
    metadata = dummy(
        input_paths=input_path,
        output_path=output_path,
        message=SECOND_MESSAGE,
        metadata=metadata,
    )
    debug(metadata)

    with open(outfile, "r") as f:
        data = json.load(f)
    assert len(data) == 2
    assert data[1]["message"] == SECOND_MESSAGE

    debug(data)
