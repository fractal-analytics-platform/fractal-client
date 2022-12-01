import asyncio

import pytest
from devtools import debug


async def test_task_collection_and_list(register_user, invoke, testdata_path):
    """
    GIVEN a pip installable package containing fractal-compatible tasks
    WHEN the collection subcommand is called
    THEN
        * the collection is initiated in the background
        * the server returns immediately and the client displays the

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
    while True:
        res1 = await invoke(f"task check-collection {state_id}")
        res1.show()
        await asyncio.sleep(1)
        if res1.data["data"]["status"] == "OK":
            break

    res2 = await invoke(f"task check-collection {state_id} --verbose")
    debug(res2)
    res2.show()
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
    while True:
        res1 = await invoke(f"task check-collection {state_id}")
        await asyncio.sleep(1)
        if res1.data["data"]["status"] == "OK":
            break

    # Second collection
    res0 = await invoke(f"task collect {PACKAGE_NAME}")
    res0.show()
    assert res0.data["data"]["info"] == "Already installed"


async def test_workflow_apply(register_user, invoke, testdata_path):
    """
    GIVEN a project and a nontrivial workflow
    WHEN the client requests to apply the workflow to the project
    THEN the workflow is scheduled and executed, and the artifacts created
    """
    PACKAGE_NAME = testdata_path / "fractal_tasks_dummy-0.1.0-py3-none-any.whl"
    WORKFLOW_NAME = "mywf"
    DATASET_NAME = "myds"

    res0 = await invoke(f"task collect {PACKAGE_NAME}")
    debug(res0)
    res0.show()
    venv_path = res0.data["data"]["venv_path"]
    debug(venv_path)
    state_id = res0.data["id"]

    while True:
        res1 = await invoke(f"task check-collection {state_id}")
        if res1.data["data"]["status"] == "OK":
            debug(res1.data)
            break
        await asyncio.sleep(1)

    res = await invoke("project new testproject prjpath")
    debug(res)
    assert res.retcode == 0

    prj = res.data
    prj_id = prj["id"]
    input_dataset_id = prj["dataset_list"][0]["id"]

    res = await invoke(f"project add-dataset {prj_id} {DATASET_NAME}")
    assert res.retcode == 0
    debug(res.data)
    output_dataset_id = res.data["id"]

    res = await invoke(
        f"dataset add-resource {prj_id} {output_dataset_id} "
        f"{testdata_path} -g 'out.json'"
    )

    res = await invoke(f"workflow new {WORKFLOW_NAME} {prj_id}")
    workflow = res.data
    workflow_id = workflow["id"]

    TASK_ID = 1
    res = await invoke(f"workflow add-task {workflow_id} {TASK_ID}")
    assert res.retcode == 0

    cmd = (
        f"workflow apply {workflow_id} {input_dataset_id} "
        f"-o {output_dataset_id} -p {prj['id']}"
    )
    debug(cmd)
    res = await invoke(cmd)
    debug(res.data)
    assert res.retcode == 0

    # TODO: verify outout


@pytest.mark.xfail
async def test_edit_task(register_user, invoke, clear_task_cache):
    # TODO:
    # Decide what it means to edit a task
    raise NotImplementedError
