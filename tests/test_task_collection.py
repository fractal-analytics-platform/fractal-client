import json
import sys
import time
from urllib.request import urlopen
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
                f"--package-extras {PACKAGE_EXTRAS} "
                "--pinned-dependency pydantic=1.10.0"
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


def test_task_collection_invalid_pinned_dependency(invoke, caplog):
    """
    Test the case where `pinned_package_versions` has the wrong format.
    """
    PACKAGE = "devtools"
    with pytest.raises(SystemExit):
        invoke(
            (
                "task collect "
                f"{PACKAGE} "
                "--pinned-dependency invalid-string"
            )
        )
    # Check that payload was prepared correctly
    log_lines = [record.message for record in caplog.records]
    debug(log_lines)
    assert "Invalid pin:" in log_lines[0]


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
    assert len(res.data) == 14


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
    debug(res0.data)
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


def test_task_collection_custom(register_user, tmp_path, invoke, caplog):
    python_interpreter = sys.executable
    package_name = "fractal-client"
    manifest = str(tmp_path / "manifest.json")

    # Download and write a valid Manifest
    manifest_url = (
        "https://github.com/fractal-analytics-platform/fractal-server/"
        "raw/main/tests/v2/fractal_tasks_mock/src/fractal_tasks_mock/"
        "__FRACTAL_MANIFEST__.json"
    )
    with urlopen(manifest_url) as f:
        manifest_dict = json.loads(f.read())
    with open(manifest, "w") as f:
        json.dump(manifest_dict, f)

    cmd = (
        f"task collect-custom --package-name {package_name} "
        f"source {python_interpreter} {manifest}"
    )
    res = invoke(cmd)
    assert res.retcode == 0
    assert isinstance(res.data, list)

    # Second API call fails (tasks with the same sources already exist)
    caplog.clear()
    with pytest.raises(SystemExit):
        res = invoke(cmd)
    # Manifest was redacted, when logging the payload
    assert '"manifest": "[value too long - redacted]"' in caplog.text

    # Missing manifest file
    cmd = (
        f"task collect-custom --package-name {package_name} "
        f"source {python_interpreter} /foo/bar"
    )
    res = invoke(cmd)
    assert res.retcode == 1
    assert "file must be on the same machine" in res.data

    cmd = (
        "--batch task collect-custom --package-root /tmp --version 2 "
        f"source2 {python_interpreter} {manifest}"
    )
    res = invoke(cmd)
    assert res.retcode == 0
    assert isinstance(res.data, str)

    # test that '--package-root' and '--package-name' are mutually exclusive
    cmd = (
        "task collect-custom "
        f"--package-root /tmp --package-name {package_name} "
        f"source3 {python_interpreter} {manifest}"
    )
    with pytest.raises(SystemExit):
        res = invoke(cmd)
