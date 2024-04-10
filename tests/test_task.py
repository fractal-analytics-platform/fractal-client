import json
import time
from urllib.request import urlretrieve

import pytest
from devtools import debug


COLLECTION_TIMEOUT = 15.0

PACKAGE_URL = (
    "https://github.com/fractal-analytics-platform/fractal-server/"
    "raw/main/tests/v2/fractal_tasks_mock/dist/"
    "fractal_tasks_mock-0.0.1-py3-none-any.whl"
)
PACKAGE_PATH = "/tmp/fractal_tasks_mock-0.0.1-py3-none-any.whl"
urlretrieve(PACKAGE_URL, PACKAGE_PATH)


def test_task_collection_command(register_user, invoke, caplog):
    """
    Test that all `task collect` options are correctly parsed and included in
    the the payload for the API request.
    """
    PACKAGE = "devtools"
    PACKAGE_VERSION = "0.11.0"
    PYTHON_VERSION = "1.2"
    PACKAGE_EXTRAS = "a,b,c"
    with pytest.raises(SystemExit):
        invoke(
            (
                "task collect "
                f"{PACKAGE} "
                f"--package-version {PACKAGE_VERSION} "
                f"--python-version {PYTHON_VERSION} "
                f"--package-extras {PACKAGE_EXTRAS}"
            )
        )

    # Check that payload was prepared correctly
    log_lines = [record.message for record in caplog.records]
    debug(log_lines)
    payload_line = next(
        line for line in log_lines if line.startswith("Original payload: ")
    )
    assert payload_line
    payload = json.loads(payload_line.strip("Original payload: "))
    debug(payload)
    assert payload["package"] == PACKAGE
    assert payload["package_version"] == PACKAGE_VERSION
    assert payload["package_extras"] == PACKAGE_EXTRAS
    assert payload["python_version"] == PYTHON_VERSION


def test_task_collection(register_user, invoke, testdata_path):
    """
    GIVEN a pip installable package containing fractal-compatible tasks
    WHEN the collection subcommand is called
    THEN
        * the collection is initiated in the background
        * the server returns immediately
    """

    res0 = invoke(f"task collect {PACKAGE_PATH}")
    debug(res0)
    res0.show()

    venv_path = res0.data["data"]["venv_path"]
    debug(venv_path)
    state_id = res0.data["id"]
    debug(state_id)

    time.sleep(0.5)

    # Wait until collection is complete
    starting_time = time.perf_counter()
    while True:
        res1 = invoke(f"task check-collection {state_id}")
        debug(res1.data)
        assert res1.retcode == 0
        res1.show()
        time.sleep(1)
        if res1.data["data"]["status"] == "OK":
            break
        assert time.perf_counter() - starting_time < COLLECTION_TIMEOUT

    res2 = invoke(f"task check-collection {state_id}" " --include-logs")
    debug(res2.data)
    assert res2.retcode == 0j
    res2.show()
    assert res2.data["data"]["status"] == "OK"

    res = invoke("task list")
    assert len(res.data) == 13


def test_repeated_task_collection(register_user, invoke, testdata_path):
    """
    GIVEN
        * a pip installable package containing fractal-compatible tasks
        * a successful collection subcommand was executed
    WHEN the collection subcommand is called a second time
    THEN
        * TBD..
    """

    res0 = invoke(f"--batch task collect {PACKAGE_PATH}")
    debug(res0)

    state_id = res0.data[0]  # extract id from batch string
    assert res0.data == "1 .fractal/fractal-tasks-mock0.0.1"

    time.sleep(0.5)

    # Wait until collection is complete
    starting_time = time.perf_counter()
    while True:
        res1 = invoke(f"task check-collection {state_id}")
        time.sleep(1)
        if res1.data["data"]["status"] == "OK":
            break
        assert time.perf_counter() - starting_time < COLLECTION_TIMEOUT

    # Second collection
    res0 = invoke(f"task collect {PACKAGE_PATH}")
    debug(res0.data)
    assert res0.data["data"]["info"] == "Already installed"


def test_task_new(register_user, invoke, tmp_path):

    # create a new task with just positional required args
    args_path = str(tmp_path / "args.json")
    args = {"image_dir": "/asdasd"}
    with open(args_path, "w") as f:
        json.dump(args, f)

    meta_path = str(tmp_path / "meta.json")
    meta = {"a": "b"}
    with open(meta_path, "w") as f:
        json.dump(meta, f)

    res = invoke(
        "task new _name  _source --command-parallel _command "
        f"--version _version --meta-parallel {meta_path} "
        f"--args-schema-parallel {args_path} "
        f"--args-schema-version 1.0.0"
    )
    debug(res.data)
    assert res.retcode == 0
    assert res.data["name"] == "_name"
    assert res.data["command_parallel"] == "_command"
    assert res.data["source"] == f"{register_user['username']}:_source"
    assert res.data["version"] == "_version"
    assert res.data["meta_parallel"] == meta
    assert res.data["args_schema_version"] == "1.0.0"

    assert "owner" in res.data.keys()
    first_task_id = int(res.data["id"])

    # create a new task with batch option
    res = invoke(
        "--batch task new _name2 _source2 --command-parallel _command2"
    )
    res.show()
    assert res.retcode == 0
    assert res.data == str(first_task_id + 1)

    # create a new task with same source as before. Note that in check_response
    # we have sys.exit(1) when status code is not the expecte one
    with pytest.raises(SystemExit) as e:
        invoke("task new _name2 _source --command-parallel _command2")
    assert e.value.code == 1

    # create a new task passing not existing file
    res = invoke(
        (
            "task new _name _source --command-parallel _command "
            "--meta-parallel ./foo.pdf"
        )
    )
    debug(res.data)
    assert res.retcode == 1
