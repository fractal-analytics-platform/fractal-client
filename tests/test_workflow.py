import json
import logging
import time
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

    # Add a WorkflowTask with --args-file and --meta-file arguments
    ARGS = {"arg": "arg_value"}
    META = {"executor": "some-executor"}
    args_file = tmp_path / "args_file.json"
    with args_file.open("w") as f:
        json.dump(ARGS, f)
    meta_file = tmp_path / "meta_file.json"
    with meta_file.open("w") as f:
        json.dump(META, f)

    cmd = f"workflow add-task {project_id} {wf.id}"
    # Test fail with no task_id nor task_name
    with pytest.raises(SystemExit):
        invoke(cmd)
    # Test fail with both task_id and task_name
    with pytest.raises(SystemExit):
        invoke(f"{cmd} --task-id {t.id} --task-name {t.name}")
    # Test fail with both task_id and version
    with pytest.raises(SystemExit):
        invoke(f"{cmd} --task-id {t.id} --task-version 1.2.3.4.5.6")
    assert caplog.records[-1].msg == (
        "Too many arguments: cannot provide both `task_id` and `task_version`."
    )

    cmd = (
        f"{cmd} --task-id {t.id} "
        f"--args-file {args_file} --meta-file {meta_file}"
    )
    debug(cmd)
    # Test success
    res = invoke(cmd)
    assert res.retcode == 0

    workflow_task = res.data
    workflow_task_id_1 = workflow_task["id"]
    debug(workflow_task)
    assert workflow_task["args"] == ARGS
    assert workflow_task["meta"] == META

    # Add a WorkflowTask by Task.name with the --batch option
    cmd = (
        f"--batch workflow add-task {project_id} {wf.id} "
        f"--task-name {t.name} --order 1"
    )
    debug(cmd)
    res = invoke(cmd)
    assert res.retcode == 0
    debug(res.data)
    workflow_task_id_2 = int(res.data)

    # Check that the WorkflowTask's in Workflow.task_list have the correct IDs
    cmd = f"workflow show {project_id} {wf.id}"
    res = invoke(cmd)
    assert res.retcode == 0
    workflow = res.data
    debug(workflow)
    list_IDs = [wftask["id"] for wftask in workflow["task_list"]]
    assert list_IDs == [workflow_task_id_1, workflow_task_id_2]


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

    cmd = f"workflow add-task {project_id} {wf.id} --task-name {task.name}"
    debug(cmd)
    res = invoke(cmd)
    assert res.retcode == 0
    debug(res.data)
    assert res.data["task"]["id"] == task.id

    # Fail when adding task via a wrong name
    with pytest.raises(SystemExit):
        cmd = (
            f"workflow add-task {project_id} {wf.id} --task-name INVALID_NAME"
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
    cmd = f"workflow add-task {project_id} {wf.id} --task-name {task1.name}"
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

    # Add task to workflow, twice
    cmd = f"workflow add-task {project_id} {wf.id} --task-id {t.id}"
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
    t = task_factory()

    # Create task, without overriding arguments
    cmd = f"workflow add-task {project_id} {wf.id} --task-id {t.id}"
    res = invoke(cmd)
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
    workflow_task_id = res.data["id"]
    cmd = (
        f"workflow edit-task {project_id} {wf.id} {workflow_task_id} "
        f"--args-file {args_file} --meta-file {meta_file}"
    )
    debug(cmd)
    res = invoke(cmd)
    assert res.retcode == 0
    assert res.data["args"] == ARGS
    assert res.data["meta"] == META

    # Check that also the workflow in the db was correctly updated
    res = invoke(f"workflow show {project_id} {wf.id}")
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
        f"workflow edit-task {project_id} {wf.id} {workflow_task_id} "
        f"--meta-file {meta_file}"
    )
    debug(cmd)
    with pytest.raises(SystemExit):
        res = invoke(cmd)


def test_workflow_apply(
    register_user, invoke, testdata_path: Path, tmp_path: Path
):
    """
    GIVEN a project and a nontrivial workflow
    WHEN the client requests to apply the workflow to the project
    THEN the workflow is scheduled and executed, and the artifacts created
    """

    # Collect tasks
    PACKAGE_NAME = (
        testdata_path / "fractal_tasks_core_alpha-0.0.1a0-py3-none-any.whl"
    )
    WORKFLOW_NAME = "mywf"
    DATASET_NAME = "myds"
    res0 = invoke(f"task collect {PACKAGE_NAME}")
    debug(res0)
    venv_path = res0.data["data"]["venv_path"]
    debug(venv_path)
    state_id = res0.data["id"]

    # Wait for task collection to end
    starting_time = time.perf_counter()
    while True:
        res1 = invoke(f"task check-collection {state_id}")
        if res1.data["data"]["status"] == "OK":
            debug(res1.data)
            break
        time.sleep(1)
        assert time.perf_counter() - starting_time < TIMEOUT

    # Create a project
    res = invoke("project new testproject")
    assert res.retcode == 0
    prj = res.data
    prj_id = prj["id"]
    res = invoke(f"project add-dataset {prj_id} test_name")
    input_dataset_id = res.data["id"]

    # Create output dataset
    res = invoke(f"project add-dataset {prj_id} {DATASET_NAME}")
    assert res.retcode == 0
    dataset = res.data
    output_dataset_id = dataset["id"]

    # Add resources to datasets
    for dataset_id in (input_dataset_id, output_dataset_id):
        res = invoke(
            f"dataset add-resource {prj_id} {dataset_id} {str(tmp_path)}"
        )
        assert res.retcode == 0

    # Create workflow and add task twice
    res = invoke(f"workflow new {WORKFLOW_NAME} {prj_id}")
    workflow = res.data
    workflow_id = workflow["id"]
    debug(workflow)
    assert res.retcode == 0
    for _ in [0, 1]:
        TASK_ID = 1
        res = invoke(
            f"workflow add-task {prj_id} {workflow_id} --task-id {TASK_ID}"
        )
        workflow_task = res.data
        debug(workflow_task)
        assert res.retcode == 0
        TASK_NAME = res.data["task"]["name"]
        debug(TASK_NAME)

    # Call `workflow apply`
    FIRST_TASK_INDEX = 0
    LAST_TASK_INDEX = 0
    cmd = (
        f"workflow apply "
        f"{prj_id} {workflow_id} {input_dataset_id} {output_dataset_id} "
        f"--start {FIRST_TASK_INDEX} --end {LAST_TASK_INDEX} "
        f'--worker-init "export SOMEVARIABLE=1"'
    )
    debug(cmd)
    res = invoke(cmd)
    job = res.data
    debug(job)
    assert res.retcode == 0
    job_id = job["id"]
    assert job["status"] == "submitted"

    # Avoid immediately calling `job show` right after `workflow apply`
    time.sleep(1)

    # Check that job completed successfully
    cmd = f"job show {prj_id} {job_id} --do-not-separate-logs"
    starting_time = time.perf_counter()
    debug(cmd)
    while True:
        res = invoke(cmd)
        job = res.data
        debug(job)
        assert res.retcode == 0
        if job["status"] == "done":
            break
        elif job["status"] == "failed":
            raise RuntimeError(job)
        time.sleep(1)
        assert time.perf_counter() - starting_time < TIMEOUT

    # Check that job has correct first_task_index and last_task_index
    # attributes
    assert job["first_task_index"] == FIRST_TASK_INDEX
    assert job["last_task_index"] == LAST_TASK_INDEX

    # Prepare and run a workflow with a failing task
    args_file = str(tmp_path / "args.json")
    with open(args_file, "w") as f:
        json.dump({"raise_error": True}, f)
    res = invoke(
        f"workflow add-task {prj_id} {workflow_id} --task-id {TASK_ID}"
        f" --args-file {args_file}"
    )
    assert res.retcode == 0
    cmd = (
        f"workflow apply "
        f"{prj_id} {workflow_id} {input_dataset_id} {output_dataset_id}"
    )
    debug(cmd)
    res = invoke(cmd)
    assert res.retcode == 0
    job_id = res.data["id"]

    # Avoid immediately calling `job show` right after `workflow apply`
    time.sleep(1)

    # Verify that status is failed, and that there is a log
    cmd = f"job show {prj_id} {job_id} --do-not-separate-logs"
    starting_time = time.perf_counter()
    while True:
        res = invoke(cmd)
        job = res.data
        debug(job)
        assert res.retcode == 0
        if job["status"] == "failed":
            break
        time.sleep(1)
        assert time.perf_counter() - starting_time < TIMEOUT
    assert job["log"] is not None

    # Prepare and submit a workflow with --batch
    res = invoke(f"workflow new OneMoreWorkflow {prj_id}")
    workflow_id = res.data["id"]
    res = invoke(
        f"workflow add-task {prj_id} {workflow_id} --task-id {TASK_ID}"
    )
    assert res.retcode == 0
    cmd = (
        f"--batch workflow apply "
        f"{prj_id} {workflow_id} {input_dataset_id} {output_dataset_id}"
    )
    debug(cmd)
    res = invoke(cmd)
    assert res.retcode == 0
    debug(res.data)
    assert isinstance(res.data, int)


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


def test_workflow_import(
    register_user,
    invoke,
    testdata_path: Path,
    task_factory,
    caplog,
):
    # collect tasks
    PACKAGE_NAME = testdata_path / "fractal_tasks_dummy-0.1.0-py3-none-any.whl"
    res0 = invoke(f"task collect {PACKAGE_NAME}")
    assert res0.retcode == 0
    state_id = res0.data["id"]
    starting_time = time.perf_counter()
    while True:
        res1 = invoke(f"task check-collection {state_id}")
        if res1.data["data"]["status"] == "OK":
            break
        time.sleep(1)
        assert time.perf_counter() - starting_time < TIMEOUT

    # create project
    PROJECT_NAME = "project_name"
    res_pj = invoke(f"project new {PROJECT_NAME}")
    assert res_pj.retcode == 0
    project_id = res_pj.data["id"]

    task_factory(source="custom_source", owner="exact-lab")

    # import workflow into project
    filename = str(testdata_path / "import-export/workflow.json")
    res = invoke(
        f"workflow import --project-id {project_id} --json-file {filename}"
    )
    debug(res.data)
    assert res.retcode == 0
    assert caplog.records[-1].msg == (
        "This workflow includes custom tasks (the ones with sources: "
        "'custom_source'), which are not meant to be portable; "
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
    debug(res.data)
    for _id in res.data.split():
        assert _id.isdigit()
    assert res.retcode == 0
    assert caplog.records[-1].msg == (
        "This workflow includes custom tasks (the ones with sources: "
        "'custom_source'), which are not meant to be portable; "
        "importing this workflow may not work as expected."
    )
