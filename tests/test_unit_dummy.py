import json
from pathlib import Path

from devtools import debug

from fractal.tasks.dummy import dummy


def test_dummy_task(tmp_path):
    input_path = [Path("/fake/input/path")]
    output_path = tmp_path

    FIRST_MESSAGE = "first run"
    outfile = dummy(
        input_paths=input_path,
        output_path=output_path,
        message=FIRST_MESSAGE,
    )

    with open(outfile, "r") as f:
        debug(f.read())
        f.seek(0)
        data = json.load(f)
    assert len(data) == 1
    assert data[0]["message"] == FIRST_MESSAGE

    SECOND_MESSAGE = "second run"
    outfile = dummy(
        input_paths=input_path,
        output_path=output_path,
        message=SECOND_MESSAGE,
    )

    with open(outfile, "r") as f:
        debug(f.read())
        f.seek(0)
        data = json.load(f)
    assert len(data) == 2
    assert data[1]["message"] == SECOND_MESSAGE

    debug(data)
