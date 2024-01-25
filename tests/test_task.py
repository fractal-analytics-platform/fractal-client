import json
import time
from pathlib import Path

import pytest
from devtools import debug

from fractal_client.cmd._aux_task_caching import TASKS_CACHE_FILENAME
from fractal_client.config import settings


COLLECTION_TIMEOUT = 15.0


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
    PACKAGE = testdata_path / "fractal_tasks_dummy-0.1.0-py3-none-any.whl"

    res0 = invoke(f"task collect {PACKAGE}")
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
        assert res1.retcode == 0
        res1.show()
        time.sleep(1)
        if res1.data["data"]["status"] == "OK":
            break
        assert time.perf_counter() - starting_time < COLLECTION_TIMEOUT

    # Add --include-logs and --do-not-separate-logs flags
    res2 = invoke(
        f"task check-collection {state_id}"
        " --include-logs"
        " --do-not-separate-logs"
    )
    assert res2.retcode == 0
    res2.show()
    assert res2.data["data"]["log"]
    assert res2.data["data"]["status"] == "OK"

    # Show again the check-collection output, without --do-not-separate-logs,
    # for visual inspection
    res = invoke(f"task check-collection {state_id} --include-logs")
    res.show()


def test_repeated_task_collection(register_user, invoke, testdata_path):
    """
    GIVEN
        * a pip installable package containing fractal-compatible tasks
        * a successful collection subcommand was executed
    WHEN the collection subcommand is called a second time
    THEN
        * TBD..
    """
    PACKAGE = testdata_path / "fractal_tasks_dummy-0.1.0-py3-none-any.whl"

    res0 = invoke(f"task collect {PACKAGE}")
    debug(res0)
    res0.show()

    venv_path = res0.data["data"]["venv_path"]
    state_id = res0.data["id"]
    debug(venv_path)
    debug(state_id)

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
    res0 = invoke(f"task collect {PACKAGE}")
    res0.show()
    assert res0.data["data"]["info"] == "Already installed"


def test_task_new(register_user, invoke, tmp_path):

    # create a new task with just positional required args

    schema_path = str(tmp_path / "schema.json")
    schema = dict(key1=1, key2=[2, "3"])
    with open(schema_path, "w") as f:
        json.dump(schema, f)

    res = invoke(
        "task new _name _command _source --version _version "
        f"--args-schema {schema_path} --args-schema-version 1.0.0"
    )
    res.show()
    assert res.retcode == 0
    assert res.data["name"] == "_name"
    assert res.data["command"] == "_command"
    assert res.data["source"] == f"{register_user['username']}:_source"
    assert res.data["input_type"] == res.data["output_type"] == "Any"
    assert res.data["version"] == "_version"
    assert res.data["meta"] == {}
    assert res.data["args_schema"] == schema
    assert res.data["args_schema_version"] == "1.0.0"

    assert "owner" in res.data.keys()
    first_task_id = int(res.data["id"])

    # create a new task with batch option
    res = invoke("--batch task new _name2 _command2 _source2")
    res.show()
    assert res.retcode == 0
    assert res.data == str(first_task_id + 1)

    # create a new task with same source as before. Note that in check_response
    # we have sys.exit(1) when status code is not the expecte one
    with pytest.raises(SystemExit) as e:
        invoke("task new _name2 _command2 _source")
    assert e.value.code == 1

    # create a new task passing not existing file
    with pytest.raises(FileNotFoundError):
        invoke("task new _name _command _source --meta-file ./foo.pdf")


def test_task_edit(
    caplog,
    register_user,
    invoke,
    invoke_as_superuser,
    tmp_path,
    testdata_path: Path,
):
    schema_path = str(tmp_path / "schema.json")
    schema = dict(key1=1, key2=[2, "3"])
    with open(schema_path, "w") as f:
        json.dump(schema, f)

    task = invoke(
        "task new _name _command _source "
        f"--args-schema {schema_path} --args-schema-version 1.0.0"
    )
    task.show()
    assert task.retcode == 0
    task_id = task.data["id"]
    NEW_NAME = "1234"

    # Test that regular user is not authorized
    with pytest.raises(SystemExit):
        res = invoke(f"task edit {task_id} --new-name {NEW_NAME}")

    # Test successful edit of string attributes
    res = invoke_as_superuser(
        f"task edit --id {task_id} --new-name {NEW_NAME}"
    )
    assert res.data["name"] == NEW_NAME
    assert res.retcode == 0

    NEW_COMMAND = "run"
    res = invoke_as_superuser(
        f"task edit --id {task_id} --new-command {NEW_COMMAND}"
    )
    assert res.data["command"] == NEW_COMMAND
    assert res.retcode == 0

    # Test fail with no task_id nor task_name
    with pytest.raises(SystemExit):
        res = invoke_as_superuser("task edit")
    # Test fail with both task_id and task_name
    with pytest.raises(SystemExit):
        res = invoke_as_superuser(
            f"task edit --id {task_id} --name {task.data['name']}"
        )
    # Test fail with both task_id and task_version
    with pytest.raises(SystemExit):
        res = invoke_as_superuser(
            f"task edit --id {task_id} --version 1.2.3.4.5.6"
        )
    assert caplog.records[-1].msg == (
        "Too many arguments: cannot provide both `id` and `version`."
    )

    NEW_TYPE = "zip"
    # Test regular updates (both by id and name)
    res = invoke_as_superuser(
        f"task edit --id {task_id} --new-input-type {NEW_TYPE}"
    )
    assert res.data["input_type"] == NEW_TYPE
    assert res.retcode == 0
    res = invoke_as_superuser(
        f"task edit --name {NEW_NAME} --new-output-type {NEW_TYPE}"
    )
    assert res.data["output_type"] == NEW_TYPE
    assert res.retcode == 0

    NEW_VERSION = "3.14"
    res = invoke_as_superuser(
        f"task edit --id {task_id} --new-version {NEW_VERSION}"
    )
    assert res.data["version"] == NEW_VERSION
    assert res.retcode == 0

    # Test regular update by name, after deleting cache
    cache_dir = Path(settings.FRACTAL_CACHE_PATH)
    cache_file = cache_dir / TASKS_CACHE_FILENAME
    cache_file.unlink(missing_ok=True)
    NEW_TYPE = "something"
    res = invoke_as_superuser(
        f"task edit --name {NEW_NAME} --new-output-type {NEW_TYPE}"
    )
    assert res.data["output_type"] == NEW_TYPE
    assert res.retcode == 0

    # Test failed update by name, after deleting cache
    cache_file.unlink(missing_ok=True)
    NEW_TYPE = "something-here"
    with pytest.raises(SystemExit):
        res = invoke_as_superuser(
            f"task edit --name INVALID_NAME --new-output-type {NEW_TYPE}"
        )

    # Test regular update by name, after creating an invalid cache
    with cache_file.open("w") as f:
        json.dump([], f)
    NEW_TYPE = "something-else"
    debug(f"task edit --name {NEW_NAME} --new-output-type {NEW_TYPE}")
    res = invoke_as_superuser(
        f"task edit --name {NEW_NAME} --new-output-type {NEW_TYPE}"
    )
    assert res.data["output_type"] == NEW_TYPE
    assert res.retcode == 0

    # Test `file not found` errors
    NEW_FILE = "foo.json"
    with pytest.raises(FileNotFoundError):
        invoke_as_superuser(f"task edit --id {task_id} --meta-file {NEW_FILE}")

    # Test successful edit of dictionary attributes
    meta_file = str(testdata_path / "task_edit_json/meta.json")
    res = invoke_as_superuser(
        f"task edit --id {task_id} --meta-file {meta_file}"
    )
    debug(res)
    debug(res.data)
    assert res.retcode == 0

    # Test succesful edit of args_schema (and args_schema_version)
    new_schema_path = str(tmp_path / "new_schema.json")
    new_schema = dict(key1=159, key3=None)
    with open(new_schema_path, "w") as f:
        json.dump(new_schema, f)
    res = invoke_as_superuser(
        f"task edit --id {task_id} --new-args-schema {new_schema_path} "
        "--new-args-schema-version 1.2.3"
    )
    assert res.retcode == 0
    assert res.data["args_schema"] == new_schema
    assert res.data["args_schema_version"] == "1.2.3"


def test_task_delete(
    register_user,
    user_factory,
    invoke,
    invoke_as_superuser,
):
    """
    Test task delete
    """
    NAME = "_name"
    VERSION = "1.0.0"
    task = invoke(f"task new {NAME} _command _source --version {VERSION}")
    task.show()
    assert task.retcode == 0
    task_id = task.data["id"]

    # Test access control
    with pytest.raises(SystemExit):
        EMAIL = "someone@example.org"
        PASSWORD = "123123"
        user_factory(email=EMAIL, password=PASSWORD)
        res = invoke(f"-u {EMAIL} -p {PASSWORD} task delete --id {task_id}")
    # Test fail "id and version"
    with pytest.raises(SystemExit):
        invoke(f"task delete --id {task_id} --version {VERSION}")
    # Test fail "name and wrong version"
    with pytest.raises(SystemExit):
        invoke(f"task delete --name {NAME} --version INVALID_VERSION")

    # Test sucess
    res = invoke("task list")
    task_list = res.data
    assert len(task_list) == 1
    res = invoke(f"task delete --name {NAME} --version {VERSION}")
    assert res.retcode == 0
    res = invoke("task list")
    task_list = res.data
    assert len(task_list) == 0


def test_task_list(register_user, invoke, testdata_path):
    """
    Install tasks from a package
    Add a custom task with an owner
    List tasks in the appropriate way
    """

    for version in ["0.1.0", "0.2.0"]:
        # Task collection
        PACKAGE = (
            testdata_path / f"fractal_tasks_dummy-{version}-py3-none-any.whl"
        )

        res = invoke(f"--batch task collect {PACKAGE}")
        assert res.retcode == 0

        state_id = res.data.split(" ")[0]

        # Wait until collection is complete
        time.sleep(1)
        starting_time = time.perf_counter()
        while True:
            res = invoke(f"task check-collection {state_id}")
            assert res.retcode == 0
            if res.data["data"]["status"] == "OK":
                break
            assert time.perf_counter() - starting_time < COLLECTION_TIMEOUT
            time.sleep(1)

    # Add custom task
    custom_task_name = "custom_task_name"
    custom_task_command = "custom_task_command"
    custom_task_source = "custom_task_source"
    custom_task_version = "9.9"
    custom_task_meta = testdata_path / "task_edit_json/meta.json"
    res = invoke(
        f"task new {custom_task_name} {custom_task_command} "
        f"{custom_task_source} --meta-file {custom_task_meta} "
        f"--version {custom_task_version}"
    )
    debug(res)
    assert res.retcode == 0

    # List tasks
    res = invoke("task list")
    assert res.retcode == 0
    task_list = res.data

    # Remove some attributes, to de-clutter output
    for task in task_list:
        for key in [
            "command",
            "meta",
            "input_type",
            "output_type",
        ]:
            if key in task:
                task.pop(key)
    debug(task_list)

    # Check that tasks are sorted as expected
    assert task_list[0]["name"] == "dummy"
    assert task_list[0]["version"] == "0.1.0"
    assert task_list[0]["owner"] is None
    assert task_list[1]["name"] == "dummy"
    assert task_list[1]["version"] == "0.2.0"
    assert task_list[1]["owner"] is None
    assert task_list[2]["name"] == "dummy parallel"
    assert task_list[2]["version"] == "0.1.0"
    assert task_list[2]["owner"] is None
    assert task_list[3]["name"] == "dummy parallel"
    assert task_list[3]["version"] == "0.2.0"
    assert task_list[3]["owner"] is None
    assert task_list[4]["name"] == custom_task_name
    assert task_list[4]["version"] == custom_task_version
    assert task_list[4]["owner"] is not None


def test_pin(register_user, invoke, testdata_path, caplog):

    PACKAGE = "fractal_tasks_core_alpha-0.0.1a0-py3-none-any.whl"
    PIN = "pydantic=1.10.3"

    with pytest.raises(SystemExit):
        invoke(
            f"task collect {testdata_path / PACKAGE} "
            f"--pinned-dependency {PIN.replace('=','~')}"
        )
    assert "Pins must be written as" in caplog.records[-1].msg

    res = invoke(
        f"task collect {testdata_path / PACKAGE} --pinned-dependency {PIN}"
    )
    assert res.retcode == 0
    state_id = res.data["id"]
    starting_time = time.perf_counter()
    while True:
        res = invoke(f"task check-collection {state_id}")
        time.sleep(1)
        if res.data["data"]["status"] == "OK":
            break
        assert time.perf_counter() - starting_time < COLLECTION_TIMEOUT

    res = invoke(f"task check-collection {state_id}")
    assert res.retcode == 0
    assert (
        f"Currently installed version of {PIN.split('=')[0]}"
        in res.data["data"]["log"]
    ) and (
        f"differs from pinned version ({PIN.split('=')[1]})"
        in res.data["data"]["log"]
    )
