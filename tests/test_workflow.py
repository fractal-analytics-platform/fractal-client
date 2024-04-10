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
    t = task_factory(type="parallel")

    INPUT_FILTERS = {"attributes": {"a": 1}, "types": {"b": True}}
    ARGS_PARALLEL = {"image_dir": "/asdasd"}
    ARGS_NON_PARALLEL = {"image_dir": "/dsadsa"}
    META_PARALLEL = {"a": "b"}
    META_NON_PARALLEL = {"c": "d"}

    input_filters_file = tmp_path / "input_filters.json"
    with input_filters_file.open("w") as f:
        json.dump(INPUT_FILTERS, f)

    args_parallel_file = tmp_path / "args_parallel_file.json"
    with args_parallel_file.open("w") as f:
        json.dump(ARGS_PARALLEL, f)

    args_non_parallel_file = tmp_path / "args_non_parallel_file.json"
    with args_non_parallel_file.open("w") as f:
        json.dump(ARGS_NON_PARALLEL, f)

    meta_parallel_file = tmp_path / "meta_paral.json"
    with meta_parallel_file.open("w") as f:
        json.dump(META_PARALLEL, f)

    meta_non_parallel_file = tmp_path / "meta_non_paral.json"
    with meta_non_parallel_file.open("w") as f:
        json.dump(META_NON_PARALLEL, f)

    cmd = f"workflow add-task {project_id} {wf.id}"
    # Test fail with no task_id nor task_name
    with pytest.raises(SystemExit):
        invoke(cmd)
    # Test fail with both task_id and task_name
    with pytest.raises(SystemExit):
        invoke(
            (
                f"{cmd} --task-id {t.id} --task-name {t.name} "
                f"--args-parallel {args_parallel_file}"
            )
        )
    # Test fail with both task_id and version
    with pytest.raises(SystemExit):
        invoke(
            (
                f"{cmd} --task-id {t.id} --task-version 1.2.3.4.5.6 "
                f"--args-parallel {args_parallel_file}"
            )
        )
    assert caplog.records[-1].msg == (
        "Too many arguments: cannot provide both `task_id` and `task_version`."
    )

    cmd_args = (
        f"{cmd} --task-id {t.id} --input-filters {input_filters_file} "
        f"--args-parallel {args_parallel_file} "
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
    assert workflow_task["args_parallel"] == ARGS_PARALLEL

    # Add a WorkflowTask by Task.name with the --batch option
    cmd_batch = (
        f"--batch workflow add-task {project_id} {wf.id} "
        f"--task-name {t.name} --order 1 --args-parallel {args_parallel_file}"
    )
    debug(cmd_batch)
    res = invoke(cmd_batch)
    assert res.retcode == 0
    debug(res.data)
    workflow_task_id_2 = int(res.data)

    # Add a WorkflowTask with meta-parallel args
    cmd_meta = (
        f"{cmd} --task-id {t.id} --input-filters {input_filters_file} "
        f"--args-parallel {args_parallel_file} "
        f"--meta-parallel {meta_parallel_file} "
    )
    debug(cmd_meta)
    # Test success
    res = invoke(cmd_meta)
    debug(res.data)
    assert res.retcode == 0

    workflow_task = res.data
    workflow_task_id_3 = workflow_task["id"]
    assert workflow_task["meta_parallel"] == META_PARALLEL
    assert workflow_task["args_parallel"] == ARGS_PARALLEL

    # Add a WorkflowTask with meta-non-parallel args
    t_non_parallel = task_factory(
        type="non_parallel", source="source non_parallel"
    )

    cmd_meta = (
        f"{cmd} --task-id {t_non_parallel.id} "
        f"--input-filters {input_filters_file} "
        f"--args-non-parallel {args_non_parallel_file} "
        f"--meta-non-parallel {meta_non_parallel_file}"
    )
    debug(cmd_meta)
    # Test success
    res = invoke(cmd_meta)
    debug(res.data)
    assert res.retcode == 0

    workflow_task = res.data
    workflow_task_id_4 = workflow_task["id"]
    assert workflow_task["meta_non_parallel"] == META_NON_PARALLEL
    assert workflow_task["args_non_parallel"] == ARGS_NON_PARALLEL

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
        workflow_task_id_4,
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
    task = task_factory(type="parallel")
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
    task1 = task_factory(type="parallel")
    task2 = task_factory(type="parallel")
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
    t = task_factory(type="parallel")

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
    t = task_factory(type="parallel")

    INPUT_FILTERS = {"attributes": {"a": 1}, "types": {"b": True}}
    ARGS_PARALLEL = {"image_dir": "/asdasd"}
    ARGS_NON_PARALLEL = {"image_dir": "/dsadsa"}
    META_PARALLEL = {"a": "b"}
    META_NON_PARALLEL = {"c": "d"}

    input_filters_file = tmp_path / "input_filters.json"
    with input_filters_file.open("w") as f:
        json.dump(INPUT_FILTERS, f)

    args_parallel_file = tmp_path / "args_parallel_file.json"
    with args_parallel_file.open("w") as f:
        json.dump(ARGS_PARALLEL, f)

    args_non_parallel_file = tmp_path / "args_non_parallel_file.json"
    with args_non_parallel_file.open("w") as f:
        json.dump(ARGS_NON_PARALLEL, f)

    meta_parallel_file = tmp_path / "meta_paral.json"
    with meta_parallel_file.open("w") as f:
        json.dump(META_PARALLEL, f)

    meta_non_parallel_file = tmp_path / "meta_non_paral.json"
    with meta_non_parallel_file.open("w") as f:
        json.dump(META_NON_PARALLEL, f)

    # Create task, without overriding arguments
    cmd = (
        f"workflow add-task {project_id} {wf.id} --task-id {t.id} "
        f"--args-parallel {args_parallel_file}"
    )
    res = invoke(cmd)
    assert res.retcode == 0

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

    # Edit workflow task
    debug(res.data)
    workflow_task_id = res.data["id"]
    cmd = (
        f"workflow edit-task {project_id} {wf.id} {workflow_task_id} "
        f"--args-parallel {args_parallel_file} "
        f"--meta-parallel {meta_parallel_file}"
    )
    debug(cmd)
    res = invoke(cmd)
    assert res.retcode == 0
    assert res.data["meta_parallel"] == META_PARALLEL
    assert res.data["args_parallel"] == ARGS_PARALLEL

    # Add a WorkflowTask with meta-non-parallel args
    t_non_parallel = task_factory(
        type="non_parallel", source="source non_parallel"
    )

    cmd = (
        f"workflow add-task {project_id} {wf.id} "
        f"--task-id {t_non_parallel.id} "
        f"--args-non-parallel {args_non_parallel_file}"
    )
    res = invoke(cmd)
    assert res.retcode == 0
    workflow_task_id = res.data["id"]

    cmd = (
        f"workflow edit-task {project_id} {wf.id} {workflow_task_id} "
        f"--input-filters {input_filters_file} "
        f"--args-non-parallel {args_non_parallel_file} "
        f"--meta-non-parallel {meta_non_parallel_file}"
    )
    debug(cmd)
    # Test success
    res = invoke(cmd)
    debug(res.data)
    assert res.retcode == 0

    workflow_task = res.data
    assert workflow_task["meta_non_parallel"] == META_NON_PARALLEL
    assert workflow_task["args_non_parallel"] == ARGS_NON_PARALLEL


def test_workflow_import(
    register_user,
    invoke,
    testdata_path: Path,
    task_factory,
    caplog,
):

    # create project
    PROJECT_NAME = "project_name"
    res_pj = invoke(f"project new {PROJECT_NAME}")
    assert res_pj.retcode == 0
    project_id = res_pj.data["id"]

    task_factory(name="task", source="PKG_SOURCE:dummy2", owner="exact-lab")

    # import workflow into project
    filename = str(testdata_path / "import-export/workflow.json")
    with open(filename, "r") as f:
        debug(f.read())
    res = invoke(
        f"workflow import --project-id {project_id} --json-file {filename}"
    )
    debug(res.data)
    assert res.retcode == 0
    assert caplog.records[-1].msg == (
        "This workflow includes custom tasks (the ones with sources: "
        "'PKG_SOURCE:dummy2'), which are not meant to be portable; "
        "importing this workflow may not work as expected."
    )

    imported_workflow = res.data

    # get the workflow from the server, and check that it is the same
    workflow_id = res.data["id"]
    res = invoke(f"workflow show {project_id} {workflow_id}")
    assert res.retcode == 0
    assert res.data == imported_workflow

    # import workflow into project, with --batch
    filename = str(testdata_path / "import-export/workflow_2.json")
    res = invoke(
        f"--batch workflow import --project-id {project_id} "
        f"--json-file {filename}"
    )
    assert res.retcode == 0
    assert res.data == "2 2"
    assert caplog.records[-1].msg == (
        "This workflow includes custom tasks (the ones with sources: "
        "'PKG_SOURCE:dummy2'), which are not meant to be portable; "
        "importing this workflow may not work as expected."
    )


def test_workflow_export(
    register_user,
    invoke,
    workflow_factory,
    tmp_path: Path,
    task_factory,
    caplog,
):

    res = invoke("project new testproject")
    assert res.retcode == 0
    project_id = res.data["id"]

    NAME = "WorkFlow"
    wf = workflow_factory(project_id=project_id, name=NAME)
    prj_id = wf.project_id
    wf_id = wf.id
    filename = str(tmp_path / "exported_wf.json")

    task = task_factory(owner="exact-lab")
    res = invoke(f"workflow add-task {prj_id} {wf_id} --task-id {task.id}")
    assert res.retcode == 0

    res = invoke(f"workflow export {prj_id} {wf_id} --json-file {filename}")
    assert res.retcode == 0
    assert caplog.records[-1].msg == (
        "This workflow includes custom tasks (the ones with sources: "
        f"'{task.source}'), which are not meant to be portable; "
        "re-importing this workflow may not work as expected."
    )
    debug(res.data)
    with open(filename, "r") as f:
        exported_wf = json.load(f)
        assert exported_wf["name"] == NAME
        assert "id" not in exported_wf
        assert "project_id" not in exported_wf
        for wftask in exported_wf["task_list"]:
            assert "id" not in wftask
            assert "task_id" not in wftask
            assert "workflow_id" not in wftask
