import asyncio
import time

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
    assert res1.data["data"]["log"] is None

    # Add --include-logs flag
    res2 = await invoke(f"task check-collection {state_id} --include-logs")
    assert res2.retcode == 0
    res2.show()
    assert res2.data["data"]["log"] is not None
    assert res2.data["data"]["status"] == "OK"

    # LIST

    res = await invoke("task list")
    res.show()
    assert res.retcode == 0
    assert len(res.data) == 2


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


@pytest.mark.xfail
async def test_edit_task(register_user, invoke, clear_task_cache):
    # TODO:
    # Decide what it means to edit a task
    raise NotImplementedError


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
