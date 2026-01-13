import json
import logging
import time
from urllib.request import urlopen

import pytest
from devtools import debug

logging.getLogger("httpx").setLevel(logging.DEBUG)

NUM_MOCK_TASKS = 17


def test_task_collection_command(invoke, caplog):
    """
    Test that all `task collect` options are correctly parsed and included in
    the the payload for the API request.
    """
    INVALID_PYTHON_VERSION = "xxx"
    PACKAGE = "devtools"
    PACKAGE_VERSION = "0.11.0"
    PYTHON_VERSION = INVALID_PYTHON_VERSION
    PACKAGE_EXTRAS = "a,b,c"
    with pytest.raises(SystemExit):
        invoke(
            "task collect "
            f"{PACKAGE} "
            f"--package-version {PACKAGE_VERSION} "
            f"--python-version {PYTHON_VERSION} "
            f"--package-extras {PACKAGE_EXTRAS} "
            "--pre-pinned-dependency pydantic=1.10.0 "
            "--post-pinned-dependency pydantic=1.10.0"
        )
    debug(caplog.text)
    assert "Server returned 422" in caplog.text
    assert f"input_value='{INVALID_PYTHON_VERSION}'" in caplog.text


def test_task_collection_invalid_pinned_dependency(invoke, caplog):
    """
    Test the case where `pre_pinned_package_versions` has the wrong format.
    """
    PACKAGE = "devtools"
    with pytest.raises(SystemExit):
        invoke(
            f"task collect {PACKAGE} --pre-pinned-dependency invalid-string"
        )
    # Check that payload was prepared correctly
    error_line = next(
        record.message
        for record in caplog.records
        if "Invalid pin:" in record.message
    )
    debug(error_line)
    assert error_line is not None


def test_task_collection(
    invoke_as_custom_user, user_factory, new_name, fractal_tasks_mock
):
    """
    GIVEN a pip installable package containing fractal-compatible tasks
    WHEN the collection subcommand is called
    THEN
        * the collection is initiated in the background
        * the server returns immediately
    """
    COLLECTION_TIMEOUT = 15.0

    new_user = dict(email=f"{new_name()}@example.org", password="1234")
    user_factory(**new_user)

    res = invoke_as_custom_user("task list", **new_user)
    initial_task_list = len(res.data)

    res0 = invoke_as_custom_user(
        f"task collect --private {fractal_tasks_mock}",
        **new_user,
    )
    debug(res0.data)
    activity_id = res0.data["id"]

    # Wait until collection is complete
    starting_time = time.perf_counter()
    while True:
        res1 = invoke_as_custom_user(
            f"task check-collection {activity_id}", **new_user
        )
        assert res1.retcode == 0
        time.sleep(0.1)
        if res1.data["status"] == "OK":
            debug(res1.data)
            break
        assert time.perf_counter() - starting_time < COLLECTION_TIMEOUT

    # Check successful status and no logs
    res2 = invoke_as_custom_user(
        f"task check-collection {activity_id}", **new_user
    )
    assert res2.retcode == 0
    assert res2.data["status"] == "OK"
    assert res2.data["log"] is None

    # Check logs
    res3 = invoke_as_custom_user(
        f"task check-collection {activity_id} --include-logs", **new_user
    )
    assert res3.retcode == 0
    assert res3.data["status"] == "OK"
    assert res3.data["log"] is not None

    # Check task list
    res = invoke_as_custom_user("task list", **new_user)
    assert len(res.data) == initial_task_list + NUM_MOCK_TASKS

    # Second collection
    with pytest.raises(SystemExit):
        invoke_as_custom_user(f"task collect {fractal_tasks_mock}", **new_user)


def test_task_collection_custom(
    user_factory, new_name, tmp_path, invoke_as_custom_user, caplog
):
    new_user = dict(email=f"{new_name()}@example.org", password="1234")
    user_factory(**new_user)

    python_interpreter = "/usr/local/bin/python"
    package_name = "pip"
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
        f"task collect-custom --private --package-name {package_name} "
        f"label {python_interpreter} {manifest}"
    )
    res = invoke_as_custom_user(cmd, **new_user)
    debug(res.data)
    assert res.retcode == 0
    assert isinstance(res.data, list)

    # Second API call fails (tasks with the same identity already exist)
    caplog.clear()
    with pytest.raises(SystemExit):
        res = invoke_as_custom_user(cmd, **new_user)
    # Manifest was redacted, when logging the payload
    assert '"manifest": "[value too long - redacted]"' in caplog.text

    # Missing manifest file
    cmd = (
        f"task collect-custom --package-name {package_name} "
        f"label {python_interpreter} /foo/bar"
    )
    res = invoke_as_custom_user(cmd, **new_user)
    assert res.retcode == 1
    assert "file must be on the same machine" in res.data

    cmd = (
        "--batch task collect-custom --private --package-root /tmp --version 2"
        f" label2 {python_interpreter} {manifest}"
    )
    res = invoke_as_custom_user(cmd, **new_user)
    assert res.retcode == 0
    assert isinstance(res.data, str)

    # test that '--package-root' and '--package-name' are mutually exclusive
    cmd = (
        "task collect-custom --private"
        f"--package-root /tmp --package-name {package_name} "
        f"label3 {python_interpreter} {manifest}"
    )
    with pytest.raises(SystemExit):
        res = invoke_as_custom_user(cmd, **new_user)
