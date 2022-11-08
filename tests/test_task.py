import asyncio

import pytest
from devtools import debug


async def test_task_collection_and_list(
    clear_db, testserver, register_user, invoke, testdata_path, temp_data_dir
):
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
    venv_path = res0.data["venv_path"]

    # cf. fixtures_testserver.py::testserver
    FRACTAL_ROOT = temp_data_dir
    collection_file = FRACTAL_ROOT / venv_path / "collection.json"

    while True:
        res1 = await invoke(f"task check-collection {venv_path}")
        debug(res1)
        res1.show()
        expected_status = "OK" if collection_file.exists() else "pending"
        assert res1.data["status"] == expected_status
        await asyncio.sleep(1)
        if expected_status == "OK":
            break

    res2 = await invoke(f"task check-collection {venv_path}")
    debug(res2)
    res2.show()
    assert res2.data["status"] == "OK"

    # LIST

    res = await invoke("task list")
    res.show()
    assert res.retcode == 0
    assert len(res.data) == 2


async def test_task_apply(
    clear_db, testserver, register_user, invoke, testdata_path
):
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
    venv_path = res0.data["venv_path"]

    while True:
        res1 = await invoke(f"task check-collection {venv_path}")
        if res1.data["status"] == "OK":
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
    raise NotImplementedError


@pytest.mark.xfail
async def test_edit_task(
    clear_db, testserver, register_user, invoke, clear_task_cache
):
    # TODO:
    # Decide what it means to edit a task
    raise NotImplementedError
