import json
import logging
from pathlib import Path

import pytest  # noqa F401
from devtools import debug

TIMEOUT = 15.0


def test_workflow_new(register_user, invoke):
    PROJECT_NAME = "project_name"
    res = invoke(f"project new {PROJECT_NAME}")
    proj = res.data
    assert proj["name"] == PROJECT_NAME
    project_id = proj["id"]

    WORKFLOW_NAME = "mywf"
    res = invoke(f"workflow new {WORKFLOW_NAME} {project_id}")
    wf = res.data
    debug(wf)
    assert res.retcode == 0
    assert wf["name"] == WORKFLOW_NAME
    assert wf["project_id"] == project_id

    # Include --batch
    WORKFLOW_NAME = "mywf-2"
    res = invoke(f"--batch workflow new {WORKFLOW_NAME} {project_id}")
    assert res.retcode == 0
    debug(res.data)
    assert isinstance(res.data, int)


def test_workflow_delete(register_user, invoke):
    # Create project
    res_pj = invoke("project new project_name")
    assert res_pj.retcode == 0
    project_id = res_pj.data["id"]

    # Create workflow
    res_wf = invoke(f"workflow new MyWorkflow {project_id}")
    workflow_id = res_wf.data["id"]
    assert res_wf.retcode == 0

    # List workflows
    res_list = invoke(f"workflow list {project_id}")
    assert res_list.retcode == 0
    debug(res_list.data)
    assert len(res_list.data) == 1

    # Delete workflow
    res_delete = invoke(f"workflow delete {project_id} {workflow_id}")
    assert res_delete.retcode == 0

    # List workflows
    res_list = invoke(f"workflow list {project_id}")
    assert res_list.retcode == 0
    debug(res_list.data)
    assert len(res_list.data) == 0


def test_workflow_edit(register_user, invoke):
    # Create a project
    res_pj = invoke("project new project_name_1")
    assert res_pj.retcode == 0
    project_id = res_pj.data["id"]

    # Add a workflow
    res_wf = invoke(f"workflow new MyWorkflow {project_id}")
    workflow_id = res_wf.data["id"]
    assert res_wf.retcode == 0

    # Fail editing with no edits
    cmd = f"workflow edit {project_id} {workflow_id}"
    with pytest.raises(SystemExit):
        res = invoke(cmd)

    # Edit workflow name
    NAME = "new-workflow-name"
    cmd = f"workflow edit {project_id} {workflow_id} --new-name {NAME}"
    debug(cmd)
    res_edit = invoke(cmd)
    assert res_edit.retcode == 0
    debug(res_edit.data)

    # List workflows, and check edit
    res = invoke(f"workflow show {project_id} {workflow_id}")
    debug(res.data)
    assert res.retcode == 0
    assert res.data["name"] == NAME


def test_workflow_list(register_user, invoke):
    PROJECT_NAME = "project_name"
    res_pj = invoke(f"project new {PROJECT_NAME}")
    project_id = res_pj.data["id"]
    debug(project_id)

    res_wf = invoke(f"workflow new WF1 {project_id}")
    res_wf.show()
    assert res_wf.retcode == 0

    res_wf = invoke(f"workflow new WF2 {project_id}")
    res_wf.show()
    assert res_wf.retcode == 0

    res_list = invoke(f"workflow list {project_id}")
    debug(res_list)
    debug(res_list.data)
    assert res_list.retcode == 0
    assert len(res_list.data) == 2


def test_workflow_list_when_two_projects_exist(register_user, invoke):
    res_pj1 = invoke("project new PRJ1")
    res_pj2 = invoke("project new PRJ2")
    project_id_1 = res_pj1.data["id"]
    project_id_2 = res_pj2.data["id"]

    NUM_WF_PROJECT_1 = 2
    NUM_WF_PROJECT_2 = 4

    for wf in range(NUM_WF_PROJECT_1):
        res_wf = invoke(f"workflow new WF{wf} {project_id_1}")
        assert res_wf.retcode == 0

    for wf in range(NUM_WF_PROJECT_2):
        res_wf = invoke(f"workflow new WF{wf} {project_id_2}")
        assert res_wf.retcode == 0

    res_list_1 = invoke(f"workflow list {project_id_1}")
    assert res_list_1.retcode == 0
    assert len(res_list_1.data) == NUM_WF_PROJECT_1

    res_list_2 = invoke(f"workflow list {project_id_2}")
    assert res_list_2.retcode == 0
    assert len(res_list_2.data) == NUM_WF_PROJECT_2


def test_workflow_add_task(
    caplog,
    invoke,
    register_user,
    task_factory,
    workflow_factory,
    tmp_path: Path,
):
    """
    GIVEN a workflow
    WHEN
        the client is invoked to add a task, with several different options
        (custom args, custom meta, --batch)
    THEN
        the WorkflowTask's are correctly registered in the db, and the returned
        object has the right properties
    """
    res = invoke("project new MyProject")
    project_id = res.data["id"]
    wf = workflow_factory(project_id=project_id)
    t = task_factory()

    INPUT_FILTERS = {"attributes": {"a": 1}, "types": {"b": True}}
    ARGS = {"image_dir": "/asdasd"}
    META = {"a": "b"}

    input_filters_file = tmp_path / "input_filters.json"
    with input_filters_file.open("w") as f:
        json.dump(INPUT_FILTERS, f)

    args_file = tmp_path / "args_file.json"
    with args_file.open("w") as f:
        json.dump(ARGS, f)

    meta_file = tmp_path / "meta.json"
    with meta_file.open("w") as f:
        json.dump(META, f)

    cmd = f"workflow add-task {project_id} {wf.id}"
    # Test fail with no task_id nor task_name
    with pytest.raises(SystemExit):
        invoke(cmd)
    # Test fail with both task_id and task_name
    with pytest.raises(SystemExit):
        invoke(
            (
                f"{cmd} --task-id {t.id} --task-name {t.name} "
                f"--args-parallel {args_file}"
            )
        )
    # Test fail with both task_id and version
    with pytest.raises(SystemExit):
        invoke(
            (
                f"{cmd} --task-id {t.id} --task-version 1.2.3.4.5.6 "
                f"--args-parallel {args_file}"
            )
        )
    assert caplog.records[-1].msg == (
        "Too many arguments: cannot provide both `task_id` and `task_version`."
    )

    cmd_args = (
        f"{cmd} --task-id {t.id} --input-filters {input_filters_file} "
        f"--args-parallel {args_file}"
    )
    debug(cmd_args)
    # Test success
    res = invoke(cmd_args)
    debug(res.data)
    assert res.retcode == 0

    workflow_task = res.data
    workflow_task_id_1 = workflow_task["id"]
    debug(workflow_task)
    assert workflow_task["input_filters"] == INPUT_FILTERS

    # Add a WorkflowTask by Task.name with the --batch option
    cmd_batch = (
        f"--batch workflow add-task {project_id} {wf.id} "
        f"--task-name {t.name} --order 1 --args-parallel {args_file}"
    )
    debug(cmd_batch)
    res = invoke(cmd_batch)
    assert res.retcode == 0
    debug(res.data)
    workflow_task_id_2 = int(res.data)

    # Add a WorkflowTask with meta-parallel args
    cmd_meta = (
        f"{cmd} --task-id {t.id} --input-filters {input_filters_file} "
        f"--args-parallel {args_file} --meta-parallel {meta_file}"
    )
    debug(cmd_meta)
    # Test success
    res = invoke(cmd_meta)
    debug(res.data)
    assert res.retcode == 0

    workflow_task = res.data
    workflow_task_id_3 = workflow_task["id"]

    # Check that the WorkflowTask's in Workflow.task_list have the correct IDs
    cmd = f"workflow show {project_id} {wf.id}"
    res = invoke(cmd)
    assert res.retcode == 0
    workflow = res.data
    debug(workflow)
    list_IDs = [wftask["id"] for wftask in workflow["task_list"]]
    assert list_IDs == [
        workflow_task_id_1,
        workflow_task_id_2,
        workflow_task_id_3,
    ]


def test_workflow_add_task_by_name(
    invoke,
    register_user,
    task_factory,
    workflow_factory,
    tmp_path: Path,
):
    """
    GIVEN a workflow and a task
    WHEN the client is invoked to add a task *by name*
    THEN the WorkflowTask is added (for a valid name) or an error is raised
    (for invalid name)
    """
    res = invoke("project new MyProject")
    project_id = res.data["id"]
    wf = workflow_factory(project_id=project_id)
    task = task_factory()
    debug(task)

    ARGS = {"image_dir": "/asdasd"}

    args_file = tmp_path / "args_file.json"
    with args_file.open("w") as f:
        json.dump(ARGS, f)

    cmd = (
        f"workflow add-task {project_id} {wf.id} --task-name {task.name} "
        f"--args-parallel {args_file}"
    )
    debug(cmd)
    res = invoke(cmd)
    assert res.retcode == 0
    debug(res.data)
    assert res.data["task"]["id"] == task.id

    # Fail when adding task via a wrong name
    with pytest.raises(SystemExit):
        cmd = (
            f"workflow add-task {project_id} {wf.id} --task-name INVALID_NAME "
            f"--args-parallel {args_file}"
        )
        debug(cmd)
        res = invoke(cmd)


@pytest.mark.skip(reason="Definition of expected behavior is ongoing")
def test_task_cache_with_non_unique_names(
    invoke,
    register_user,
    task_factory,
    workflow_factory,
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
):
    """
    GIVEN two tasks with the same name
    WHEN the client is invoked to list the tasks
    THEN
        * A warning is raised that the cache won't be written
        * Addressing tasks by name raises a FileNotFoundError
    """

    res = invoke("project new MyProject")
    project_id = res.data["id"]
    ARGS = {"image_dir": "/asdasd"}

    args_file = tmp_path / "args_file.json"
    with args_file.open("w") as f:
        json.dump(ARGS, f)
    # Create two tasks with the same name
    task1 = task_factory()
    task2 = task_factory()
    assert task1.name == task2.name

    # Verify that a warning is raised upon creating the cache file
    caplog.set_level(logging.WARNING)
    res = invoke("task list")
    assert res.retcode == 0
    debug(caplog.text)
    assert "Cannot" in caplog.text

    # Verify that adding tasks to a worfklow by name (as opposed to "by id")
    # fails because of missing cache file
    wf = workflow_factory(project_id=project_id)
    cmd = (
        f"workflow add-task {project_id} {wf.id} --task-name {task1.name} "
        f"--args-parallel {args_file}"
    )
    debug(cmd)
    with pytest.raises(FileNotFoundError):
        res = invoke(cmd)


def test_workflow_rm_task(
    invoke,
    register_user,
    task_factory,
    workflow_factory,
    tmp_path: Path,
):
    # Create project, workflow and task
    res = invoke("project new MyProject")
    project_id = res.data["id"]
    wf = workflow_factory(project_id=project_id)
    t = task_factory()

    ARGS = {"image_dir": "/asdasd"}

    args_file = tmp_path / "args_file.json"
    with args_file.open("w") as f:
        json.dump(ARGS, f)

    # Add task to workflow, twice
    cmd = (
        f"workflow add-task {project_id} {wf.id} --task-id {t.id} "
        f"--args-parallel {args_file}"
    )
    res = invoke(cmd)
    assert res.retcode == 0
    res = invoke(cmd)
    assert res.retcode == 0
    workflow_task_id_1 = res.data["id"]

    # Remove task 1 from workflow
    cmd = f"workflow rm-task {project_id} {wf.id} {workflow_task_id_1}"
    debug(cmd)
    res = invoke(cmd)
    assert res.retcode == 0
    debug(res.data)


@pytest.mark.skip(reason="Skip until server will implement this feature")
def test_workflow_edit_task(
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

    res = invoke("project new MyProject")
    project_id = res.data["id"]
    wf = workflow_factory(project_id=project_id)
    t = task_factory()
    ARGS = {"image_dir": "/asdasd"}

    args_file = tmp_path / "args_file.json"
    with args_file.open("w") as f:
        json.dump(ARGS, f)
    # Create task, without overriding arguments
    cmd = (
        f"workflow add-task {project_id} {wf.id} --task-id {t.id} "
        f"--args-parallel {args_file}"
    )
    res = invoke(cmd)
    assert res.retcode == 0

    INPUT_FILTERS = {"attributes": {"a": 1}, "types": {"b": True}}

    input_filters_file = tmp_path / "input_filters_file.json"
    with input_filters_file.open("w") as f:
        json.dump(INPUT_FILTERS, f)

    # Edit workflow task
    debug(res.data)
    workflow_task_id = res.data["id"]
    cmd = (
        f"workflow edit-task {project_id} {wf.id} {workflow_task_id} "
        f"--input-filters {input_filters_file}"
    )
    debug(cmd)
    res = invoke(cmd)
    assert res.retcode == 0
    assert res.data["input_filters"] == INPUT_FILTERS
