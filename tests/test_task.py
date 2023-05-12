import asyncio
import time
from pathlib import Path

import pytest
from devtools import debug


COLLECTION_TIMEOUT = 15.0


async def test_task_collection_and_list(register_user, invoke, testdata_path):
    """
    GIVEN a pip installable package containing fractal-compatible tasks

    WHEN the collection subcommand is called
    THEN
        * the collection is initiated in the background
        * the server returns immediately

    WHEN the list command is called
    THEN the tasks collected are shown
    """
    PACKAGE_NAME = testdata_path / "fractal_tasks_dummy-0.1.0-py3-none-any.whl"

    res0 = await invoke(f"task collect {PACKAGE_NAME}")
    debug(res0)
    res0.show()

    venv_path = res0.data["data"]["venv_path"]
    debug(venv_path)
    state_id = res0.data["id"]
    debug(state_id)

    # Wait until collection is complete
    starting_time = time.perf_counter()
    while True:
        res1 = await invoke(f"task check-collection {state_id}")
        assert res1.retcode == 0
        res1.show()
        await asyncio.sleep(1)
        if res1.data["data"]["status"] == "OK":
            break
        assert time.perf_counter() - starting_time < COLLECTION_TIMEOUT

    # Add --include-logs and --do-not-separate-logs flags
    res2 = await invoke(
        f"task check-collection {state_id}"
        " --include-logs"
        " --do-not-separate-logs"
    )
    assert res2.retcode == 0
    res2.show()
    assert res2.data["data"]["log"]
    assert res2.data["data"]["status"] == "OK"

    # List tasks
    res = await invoke("task list")
    res.show()
    assert res.retcode == 0
    assert len(res.data) == 2

    # Show again the check-collection output, without --do-not-separate-logs,
    # for visual inspection
    res = await invoke(f"task check-collection {state_id} --include-logs")
    res.show()


async def test_repeated_task_collection(register_user, invoke, testdata_path):
    """
    GIVEN
        * a pip installable package containing fractal-compatible tasks
        * a successful collection subcommand was executed
    WHEN the collection subcommand is called a second time
    THEN
        * TBD..
    """
    PACKAGE_NAME = testdata_path / "fractal_tasks_dummy-0.1.0-py3-none-any.whl"

    res0 = await invoke(f"task collect {PACKAGE_NAME}")
    debug(res0)
    res0.show()

    venv_path = res0.data["data"]["venv_path"]
    state_id = res0.data["id"]
    debug(venv_path)
    debug(state_id)

    # Wait until collection is complete
    starting_time = time.perf_counter()
    while True:
        res1 = await invoke(f"task check-collection {state_id}")
        await asyncio.sleep(1)
        if res1.data["data"]["status"] == "OK":
            break
        assert time.perf_counter() - starting_time < COLLECTION_TIMEOUT

    # Second collection
    res0 = await invoke(f"task collect {PACKAGE_NAME}")
    res0.show()
    assert res0.data["data"]["info"] == "Already installed"


async def test_task_new(register_user, invoke):

    # create a new task with just positional required args
    res = await invoke("task new _name _command _source")
    res.show()
    assert res.retcode == 0
    assert res.data["name"] == "_name"
    assert res.data["command"] == "_command"
    assert res.data["source"] == "_source"
    assert res.data["input_type"] == res.data["output_type"] == "Any"
    assert res.data["default_args"] == res.data["meta"] == {}
    first_task_id = int(res.data["id"])

    # create a new task with batch option
    res = await invoke("--batch task new _name2 _command2 _source2")
    res.show()
    assert res.retcode == 0
    assert res.data == str(first_task_id + 1)

    # create a new task with same source as before. Note that in check_response
    # we have sys.exit(1) when status code is not the expecte one
    with pytest.raises(SystemExit) as e:
        await invoke("task new _name2 _command2 _source")
    assert e.value.code == 1

    # create a new task passing not existing file
    with pytest.raises(FileNotFoundError):
        await invoke("task new _name _command _source --meta-file ./foo.pdf")


async def test_task_edit(
    register_user,
    invoke,
    invoke_as_superuser,
    testdata_path: Path,
):
    task = await invoke("task new _name _command _source")
    task.show()
    assert task.retcode == 0
    task_id = task.data["id"]
    NEW = "new"

    # Test that regular user is not authorized
    with pytest.raises(SystemExit):
        res = await invoke(f"task edit {task_id} --name {NEW}")

    # Test successful edit of string attributes
    res = await invoke_as_superuser(f"task edit {task_id} --name {NEW}")
    assert res.data["name"] == NEW
    assert res.retcode == 0
    res = await invoke_as_superuser(f"task edit {task_id} --command {NEW}")
    assert res.data["command"] == NEW
    assert res.retcode == 0
    res = await invoke_as_superuser(f"task edit {task_id} --input-type {NEW}")
    assert res.data["input_type"] == NEW
    assert res.retcode == 0
    res = await invoke_as_superuser(f"task edit {task_id} --output-type {NEW}")
    assert res.data["output_type"] == NEW
    assert res.retcode == 0

    # Test `file not found` errors
    with pytest.raises(FileNotFoundError):
        res = await invoke_as_superuser(
            f"task edit {task_id} --default-args-file {NEW}"
        )
    with pytest.raises(FileNotFoundError):
        await invoke_as_superuser(f"task edit {task_id} --meta-file {NEW}")

    # Test successful edit of dictionary attributes
    args_file = str(testdata_path / "task_edit_json/default_args.json")
    res = await invoke_as_superuser(
        f"task edit {task_id} --default-args-file {args_file}"
    )
    debug(res)
    debug(res.data)
    assert res.retcode == 0
    meta_file = str(testdata_path / "task_edit_json/meta.json")
    res = await invoke_as_superuser(
        f"task edit {task_id} --meta-file {meta_file}"
    )
    debug(res)
    debug(res.data)
    assert res.retcode == 0


async def test_task_delete(register_user, invoke):
    """
    This is currently a placeholder test, since task-delete is not implemented
    """
    res = await invoke("task new _name _command _source")
    debug(res.data)
    assert res.retcode == 0
    task_id = res.data["id"]
    with pytest.raises(NotImplementedError):
        debug(f"task delete {task_id}")
        res = await invoke(f"task delete {task_id}")
