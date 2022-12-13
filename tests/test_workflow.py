import asyncio
import json
import logging
import time
from pathlib import Path

import pytest  # noqa F401
from devtools import debug

TIMEOUT = 15.0


async def test_workflow_new(register_user, invoke):
    PROJECT_NAME = "project_name"
    PROJECT_PATH = "project_path"
    WORKFLOW_NAME = "mywf"
    res_pj = await invoke(f"project new {PROJECT_NAME} {PROJECT_PATH}")
    assert res_pj.data["name"] == PROJECT_NAME
    assert res_pj.data["project_dir"] == PROJECT_PATH

    res_wf = await invoke(f"workflow new {WORKFLOW_NAME} {res_pj.data['id']}")
    res_wf.show()
    assert res_wf.retcode == 0
    assert res_wf.data["name"] == WORKFLOW_NAME
    assert res_wf.data["project_id"] == res_pj.data["id"]


async def test_workflow_list(register_user, invoke):
    PROJECT_NAME = "project_name"
    PROJECT_PATH = "project_path"
    res_pj = await invoke(f"project new {PROJECT_NAME} {PROJECT_PATH}")
    project_id = res_pj.data["id"]
    debug(project_id)

    res_wf = await invoke(f"workflow new WF1 {project_id}")
    res_wf.show()
    assert res_wf.retcode == 0

    res_wf = await invoke(f"workflow new WF2 {project_id}")
    res_wf.show()
    assert res_wf.retcode == 0

    res_list = await invoke(f"workflow list {project_id}")
    debug(res_list)
    debug(res_list.data)
    assert res_list.retcode == 0
    assert len(res_list.data) == 2


async def test_workflow_list_when_two_projects_exist(
    register_user, invoke, tmp_path: Path
):
    res_pj1 = await invoke(f"project new PRJ1 {str(tmp_path)}/prj1")
    res_pj2 = await invoke(f"project new PRJ2 {str(tmp_path)}/prj2")
    project_id_1 = res_pj1.data["id"]
    project_id_2 = res_pj2.data["id"]

    NUM_WF_PROJECT_1 = 2
    NUM_WF_PROJECT_2 = 4

    for wf in range(NUM_WF_PROJECT_1):
        res_wf = await invoke(f"workflow new WF{wf} {project_id_1}")
        assert res_wf.retcode == 0

    for wf in range(NUM_WF_PROJECT_2):
        res_wf = await invoke(f"workflow new WF{wf} {project_id_2}")
        assert res_wf.retcode == 0

    res_list_1 = await invoke(f"workflow list {project_id_1}")
    assert res_list_1.retcode == 0
    assert len(res_list_1.data) == NUM_WF_PROJECT_1

    res_list_2 = await invoke(f"workflow list {project_id_2}")
    assert res_list_2.retcode == 0
    assert len(res_list_2.data) == NUM_WF_PROJECT_2


async def test_add_task(
    invoke,
    register_user,
    task_factory,
    workflow_factory,
    tmp_path: Path,
):
    """
    GIVEN a workflow
    WHEN the client is invoked to add a task, including custom args
    THEN
        the WorkflowTask is correctly registered in the db, including custom
        args
    """
    wf = await workflow_factory()
    t = await task_factory()

    ARGS = {"arg": "arg_value"}
    META = {"executor": "some-executor"}
    args_file = tmp_path / "args_file.json"
    with args_file.open("w") as f:
        json.dump(ARGS, f)
    meta_file = tmp_path / "meta_file.json"
    with meta_file.open("w") as f:
        json.dump(META, f)

    cmd = (
        f"workflow add-task {wf.id} {t.id} "
        f"--args-file {args_file} --meta-file {meta_file}"
    )
    debug(cmd)
    res = await invoke(cmd)
    assert res.retcode == 0
    debug(res.data)
    assert res.data["task_list"][0]["args"] == ARGS
    assert res.data["task_list"][0]["meta"] == META


async def test_add_task_by_name(
    invoke,
    register_user,
    task_factory,
    workflow_factory,
    tmp_path: Path,
    clear_task_cache,
):
    """
    GIVEN a workflow and a task
    WHEN the client is invoked to add a task *by name*
    THEN the WorkflowTask is correctly registered in the db
    """
    wf = await workflow_factory()
    task = await task_factory()
    debug(task)

    cmd = f"workflow add-task {wf.id} {task.name}"
    debug(cmd)
    res = await invoke(cmd)
    assert res.retcode == 0
    debug(res)
    debug(res.data)
    assert res.data["task_list"][0]["id"] == task.id


async def test_task_cache_fails_for_non_unique_names(
    invoke,
    register_user,
    task_factory,
    workflow_factory,
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
    clear_task_cache,
):
    """
    GIVEN two tasks with the same name
    WHEN the client is invoked to list the tasks
    THEN
        * A warning is raised that the cache won't be written
        * Addressing tasks by name raises a FileNotFoundError
    """

    # Create two tasks with the same name
    task1 = await task_factory()
    task2 = await task_factory()
    assert task1.name == task2.name

    # Verify that a warning is raised upon creating the cache file
    caplog.set_level(logging.WARNING)
    res = await invoke("task list")
    assert res.retcode == 0
    debug(caplog.text)
    assert "Cannot write task-list cache" in caplog.text

    # Verify that adding tasks to a worfklow by name (as opposed to "by id")
    # fails because of missing cache file
    wf = await workflow_factory()
    cmd = f"workflow add-task {wf.id} {task1.name}"
    debug(cmd)
    with pytest.raises(FileNotFoundError):
        res = await invoke(cmd)


async def test_edit_workflow_task(
    invoke,
    register_user,
    task_factory,
    workflow_factory,
    tmp_path: Path,
):
    """
    GIVEN a workflow
    WHEN the client is invoked to add a task, including custom args
    THEN
        the WorkflowTask is correctly registered in the db, including custom
        gargs
    """
    wf = await workflow_factory()
    t = await task_factory()

    # Create task, without overriding arguments
    cmd = f"workflow add-task {wf.id} {t.id}"
    res = await invoke(cmd)
    assert res.retcode == 0

    # New arguments to be used
    ARGS = {"some_arg": "some_value"}
    META = {"executor": "cpu-low"}

    args_file = tmp_path / "args_file.json"
    with args_file.open("w") as f:
        json.dump(ARGS, f)
    meta_file = tmp_path / "meta_file.json"
    with meta_file.open("w") as f:
        json.dump(META, f)

    # Edit workflow task
    debug(res.data)
    workflow_task_id = res.data["task_list"][0]["id"]
    cmd = (
        f"workflow edit-task {wf.id} {workflow_task_id} "
        f"--args-file {args_file} --meta-file {meta_file}"
    )
    debug(cmd)
    res = await invoke(cmd)
    assert res.retcode == 0
    assert res.data["args"] == ARGS
    assert res.data["meta"] == META

    # Check that also the workflow in the db was correctly updated
    res = await invoke(f"workflow show {wf.id}")
    assert res.retcode == 0
    debug(res.data)
    assert res.data["task_list"][0]["args"] == ARGS
    assert res.data["task_list"][0]["meta"] == META

    # Check if the correct error is raised where parallelization_level
    # is set
    META_err = {"parallelization_level": "XXX"}
    meta_file = tmp_path / "meta_file.json"
    with meta_file.open("w") as f:
        json.dump(META_err, f)

    workflow_task_id = res.data["task_list"][0]["id"]
    cmd = (
        f"workflow edit-task {wf.id} {workflow_task_id} "
        f"--meta-file {meta_file}"
    )
    debug(cmd)
    with pytest.raises(ValueError):
        res = await invoke(cmd)


async def test_workflow_apply(
    register_user, invoke, testdata_path: Path, tmp_path: Path
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
    venv_path = res0.data["data"]["venv_path"]
    debug(venv_path)
    state_id = res0.data["id"]

    starting_time = time.perf_counter()
    while True:
        res1 = await invoke(f"task check-collection {state_id}")
        if res1.data["data"]["status"] == "OK":
            debug(res1.data)
            break
        await asyncio.sleep(1)
        assert time.perf_counter() - starting_time < TIMEOUT

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
    debug(res.data)
    TASK_NAME = res.data["task_list"][0]["task"]["name"]
    debug(TASK_NAME)

    cmd = (
        f"workflow apply {workflow_id} {input_dataset_id} "
        f"-o {output_dataset_id} -p {prj['id']}"
    )
    debug(cmd)
    res = await invoke(cmd)
    debug(res.data)
    job_id = res.data["id"]
    assert res.retcode == 0
    assert res.data["status"] == "submitted"

    # Check that job completed successfully
    cmd = f"workflow job-status --job-id {job_id}"
    starting_time = time.perf_counter()
    debug(cmd)
    while True:
        res = await invoke(cmd)
        assert res.retcode == 0
        if res.data["status"] == "done":
            break
        await asyncio.sleep(1)
        assert time.perf_counter() - starting_time < TIMEOUT
    debug(res.data)
    assert res.data["history"][0] == TASK_NAME

    cmd = f"--batch workflow job-status --job-id {job_id}"
    res = await invoke(cmd)
    assert res.retcode == 0
    assert res.data == "done"

    # Prepare and run a workflow with a failing task
    args_file = str(tmp_path / "args.json")
    with open(args_file, "w") as f:
        json.dump({"raise_error": True}, f)
    res = await invoke(
        f"workflow add-task {workflow_id} {TASK_ID}"
        f" --args-file {args_file}"
    )
    assert res.retcode == 0
    cmd = (
        f"workflow apply {workflow_id} {input_dataset_id} "
        f"-o {output_dataset_id} -p {prj['id']}"
    )
    debug(cmd)
    res = await invoke(cmd)
    assert res.retcode == 0
    job_id = res.data["id"]

    # Verify that status is failed, and that there is a log
    cmd = f"workflow job-status --job-id {job_id}"
    res = await invoke(cmd)
    assert res.retcode == 0
    debug(res.data)
    assert res.data["status"] == "failed"
    assert res.data["log"] is not None
    # Note: the failing task is not added to the history
    assert len(res.data["history"]) > 0
    print(res.data["log"])
