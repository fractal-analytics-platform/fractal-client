import json
import time
from pathlib import Path
from urllib.request import urlretrieve

import pytest
from devtools import debug

from fractal_client.cmd._aux_task_caching import TASKS_CACHE_FILENAME
from fractal_client.config import settings


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


def test_task_edit(
    caplog,
    register_user,
    invoke,
    invoke_as_superuser,
    tmp_path,
    testdata_path: Path,
):

    args_path = str(tmp_path / "args.json")
    args = {"image_dir": "/asdasd"}
    with open(args_path, "w") as f:
        json.dump(args, f)

    meta_path = str(tmp_path / "meta.json")
    meta = {"a": "b"}
    with open(meta_path, "w") as f:
        json.dump(meta, f)

    task = invoke(
        "task new _name  _source --command-parallel _command "
        f"--version _version --meta-parallel {meta_path} "
        f"--args-schema-parallel {args_path} "
        f"--args-schema-version 1.0.0"
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

    NEW_COMMAND_PARALLEL = "run_parallel"
    res = invoke_as_superuser(
        (
            f"task edit --id {task_id} "
            f"--command-parallel {NEW_COMMAND_PARALLEL}"
        )
    )
    assert res.data["command_parallel"] == NEW_COMMAND_PARALLEL
    assert res.retcode == 0

    # Add non-parallel task and test command-non-parallel

    meta_path = str(tmp_path / "meta.json")
    meta = {"a": "b"}
    with open(meta_path, "w") as f:
        json.dump(meta, f)

    task_np = invoke(
        (
            f"task new _name_np _source_np "
            f"--command-non-parallel _command_np "
            f"--version 1.0.1 --meta-non-parallel {meta_path}"
        )
    )

    NEW_COMMAND_NON_PARALLEL = "run_non_parallel"
    res = invoke_as_superuser(
        (
            f"task edit --id {task_np.data['id']} "
            f"--command-non-parallel {NEW_COMMAND_NON_PARALLEL}"
        )
    )
    assert res.data["command_non_parallel"] == NEW_COMMAND_NON_PARALLEL
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
    # Test fail "name and wrong version"
    with pytest.raises(SystemExit):
        invoke("task delete --name INVALID_NAME --version INVALID_VERSION")

    input_types = {"input": True, "output": False}

    i_types_path = str(tmp_path / "itypes.json")
    with open(i_types_path, "w") as f:
        json.dump(input_types, f)

    output_types = {"input": False, "output": True}

    o_types_path = str(tmp_path / "otypes.json")
    with open(o_types_path, "w") as f:
        json.dump(output_types, f)

    # Test regular updates (both by id and name)
    res = invoke_as_superuser(
        f"task edit --id {task_id} --input-types {i_types_path}"
    )
    assert res.data["input_types"] == input_types
    assert res.retcode == 0
    res = invoke_as_superuser(
        f"task edit --name {NEW_NAME} --output-types {o_types_path}"
    )
    assert res.data["output_types"] == output_types
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

    res = invoke_as_superuser(
        f"task edit --name {NEW_NAME} --output-types {o_types_path}"
    )
    assert res.data["output_types"] == output_types
    assert res.retcode == 0

    # Test failed update by name, after deleting cache
    cache_file.unlink(missing_ok=True)

    fail_output_types = {"input": True, "output": False}

    f_o_types_path = str(tmp_path / "fotypes.json")
    with open(f_o_types_path, "w") as f:
        json.dump(fail_output_types, f)

    with pytest.raises(SystemExit):
        res = invoke_as_superuser(
            f"task edit --name INVALID_NAME --output-type {f_o_types_path}"
        )

    # Test regular update by name, after creating an invalid cache
    with cache_file.open("w") as f:
        json.dump([], f)

    new_output_types = {"input": False, "output": True}

    n_o_types_path = str(tmp_path / "notypes.json")
    with open(n_o_types_path, "w") as f:
        json.dump(new_output_types, f)

    res = invoke_as_superuser(
        f"task edit --name {NEW_NAME} --output-types {n_o_types_path}"
    )
    assert res.data["output_types"] == new_output_types
    assert res.retcode == 0


def test_task_delete(
    register_user,
    user_factory,
    invoke,
    invoke_as_superuser,
    tmp_path,
):
    """
    Test task delete
    """
    NAME = "_name"
    VERSION = "1.0.0"

    meta_path = str(tmp_path / "meta.json")
    meta = {"a": "b"}
    with open(meta_path, "w") as f:
        json.dump(meta, f)

    task = invoke(
        (
            f"task new {NAME}  _source --command-parallel _command "
            f"--version {VERSION} --meta-parallel {meta_path}"
        )
    )

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

    # Test success
    res = invoke("task list")
    task_list = res.data
    assert len(task_list) == 1
    res = invoke(f"task delete --name {NAME} --version {VERSION}")
    assert res.retcode == 0
    res = invoke("task list")
    task_list = res.data
    assert len(task_list) == 0
