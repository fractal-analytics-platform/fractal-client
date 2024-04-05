import json
import time
from pathlib import Path
from urllib.request import urlretrieve

import pytest  # noqa F401
from devtools import debug
from fractal_server.utils import get_timestamp

TIMEOUT = 15.0

LOG = "Here are some logs"


@pytest.mark.parametrize("status", ["done", "failed"])
def test_job_show(
    register_user,
    invoke,
    tmp_path: Path,
    status: str,
    workflow_factory,
    project_factory,
    job_factory,
):
    """
    GIVEN a job entry in the database
    WHEN calling the `job show` command with multiple options
    THEN the client response has the expected status and log attributes
    """

    # Create mock Workflow/ApplyWorkflow objects
    res = invoke("project new prj0")
    project_id = res.data["id"]
    wf = workflow_factory(project_id=project_id)
    debug(wf)
    log = LOG
    workflow_path = tmp_path / f"workflow_{wf.id}"
    workflow_path.mkdir()
    job = job_factory(
        working_dir=workflow_path.as_posix(),
        worfklow_id=wf.id,
        status=status,
        log=log,
        workflow_dump={
            "name": "my workflow",
            "id": 1,
            "project_id": 1,
            "task_list": [],
            "timestamp_created": str(get_timestamp()),
        },
    )
    debug(job)

    # Check `job show` output
    cmd = f"job show {project_id} {job.id}"
    debug(cmd)
    res = invoke(cmd)
    assert res.retcode == 0
    assert res.data["status"] == status
    debug(res.data)
    assert LOG in res.data["log"]
    res.show()

    # Check `job show` output with --batch
    cmd = f"--batch job show {project_id} {job.id}"
    res = invoke(cmd)
    assert res.retcode == 0
    assert res.data == status


def test_job_list(
    register_user,
    invoke,
    tmp_path: Path,
    project_factory,
    workflow_factory,
    job_factory,
):
    """
    GIVEN several job entries in the database
    WHEN calling the `job list ` command
    THEN the client response lists the jobs associated to the project
    """

    # Create mock Project/Workflow/ApplyWorkflow objects
    res = invoke("project new prj0")
    project_id = res.data["id"]
    wf_1 = workflow_factory(project_id=project_id)
    wf_2 = workflow_factory(project_id=project_id)
    wd_1 = tmp_path / f"workflow_{wf_1.id}"
    wd_2 = tmp_path / f"workflow_{wf_2.id}"
    job1 = job_factory(
        working_dir=wd_1.as_posix(),
        worfklow_id=wf_1.id,
        status="running",
        workflow_dump={
            "name": "my workflow",
            "id": 1,
            "project_id": 1,
            "task_list": [],
            "timestamp_created": str(get_timestamp()),
        },
    )
    job2 = job_factory(
        working_dir=wd_2.as_posix(),
        worfklow_id=wf_2.id,
        status="done",
        workflow_dump={
            "name": "my workflow",
            "id": 1,
            "project_id": 1,
            "task_list": [],
            "timestamp_created": str(get_timestamp()),
        },
    )
    debug(job1)
    debug(job2)

    # Check `job list` output with --batch option
    cmd = f"--batch job list {project_id}"
    debug(cmd)
    res = invoke(cmd)
    debug(res.data)
    assert res.retcode == 0
    job_ids = [int(i) for i in res.data.split()]
    assert job1.id in job_ids
    assert job2.id in job_ids

    # Check `job list` output
    cmd = f"job list {project_id}"
    debug(cmd)
    res = invoke(cmd)
    assert res.retcode == 0
    debug(res.data)
    # There is not much to assert here, apart from successful invocation of the
    # command. We add a res.show() for when pytest is run with the -s flag
    res.show()


def test_job_download_logs(
    register_user,
    invoke,
    tmp_path: Path,
    workflow_factory,
    job_factory,
):
    """
    Test the `job download-logs` command
    """

    # Create mock Workflow/ApplyWorkflow objects
    res = invoke("project new prj0")
    project_id = res.data["id"]
    wf = workflow_factory(project_id=project_id)
    wd = tmp_path / f"workflow_{wf.id}"
    job = job_factory(
        working_dir=wd.as_posix(),
        worfklow_id=wf.id,
        status="running",
    )

    # Write something in a logfile within the workflow directory
    LOGFILE = "log.txt"
    wd.mkdir()
    with (wd / LOGFILE).open("w") as f:
        f.write(LOG)

    # Check that download-logs fails if the output folder already exists
    output_fail = tmp_path / "output_dir_for_logs_fail"
    output_fail.mkdir()
    cmd = (
        f"job download-logs {job.project_id} {job.id} "
        f"--output {str(output_fail)}"
    )
    debug(cmd)
    res = invoke(cmd)
    assert res.retcode == 1

    # Check standard `job download-logs` output
    output = tmp_path / "output_dir_for_logs"
    cmd = (
        f"job download-logs {job.project_id} {job.id} "
        f"--output {str(output)}"
    )
    debug(cmd)
    res = invoke(cmd)
    assert res.retcode == 0
    debug(res.data)

    # Check that the logfile is there and has the right content
    logfile = output / LOGFILE
    assert logfile.exists()
    with logfile.open("r") as f:
        contents = f.read()
    assert contents == LOG


def test_job_stop(
    register_user,
    invoke,
    tmp_path: Path,
    workflow_factory,
    project_factory,
    job_factory,
    caplog,
):
    """
    GIVEN a job entry in the database
    WHEN calling the `job stop` command
    THEN the client response has the expected status

    NOTE 1: This is a test of the client command, not a test of the
    corresponding fractal-server feature

    NOTE 2: We don't have a fractal-server instance with SLURM backend in the
    fractal client CI, so this command can be tested for consistency but not
    for functionality.
    """

    # Create mock Workflow/ApplyWorkflow objects
    res = invoke("project new prj0")
    project_id = res.data["id"]
    wf = workflow_factory(project_id=project_id)
    workflow_path = tmp_path / f"workflow_{wf.id}"
    workflow_path.mkdir()
    job = job_factory(
        working_dir=workflow_path.as_posix(),
        worfklow_id=wf.id,
    )
    debug(job)

    # Call `job stop` (this will fail because FRACTAL_RUNNER_BACKEND="local"
    cmd = f"job stop {project_id} {job.id}"
    debug(cmd)
    with pytest.raises(SystemExit):
        res = invoke(cmd)
    debug(caplog.text)
    assert "Stopping a job execution is not implemented" in caplog.text
    caplog.clear()


def test_job_submit(
    register_user, invoke, testdata_path: Path, tmp_path: Path
):
    """
    GIVEN a project and a nontrivial workflow
    WHEN the client requests to apply the workflow to the project
    THEN the workflow is scheduled and executed, and the artifacts created
    """

    # Collect tasks
    PACKAGE_URL = (
        "https://github.com/fractal-analytics-platform/fractal-server/"
        "raw/main/tests/v2/fractal_tasks_mock/dist/"
        "fractal_tasks_mock-0.0.1-py3-none-any.whl"
    )
    PACKAGE_PATH = "/tmp/fractal_tasks_mock-0.0.1-py3-none-any.whl"
    urlretrieve(PACKAGE_URL, PACKAGE_PATH)

    WORKFLOW_NAME = "mywf"
    res0 = invoke(f"task collect {PACKAGE_PATH}")
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
    res = invoke(f"project add-dataset {prj_id} test_name /tmp")
    dataset_id = res.data["id"]

    # Create workflow and add task twice
    res = invoke(f"workflow new {WORKFLOW_NAME} {prj_id}")
    workflow = res.data
    workflow_id = workflow["id"]
    args_file = str(tmp_path / "args_file.json")
    with open(args_file, "w") as f:
        json.dump({"image_dir": "/asdasd"}, f)
    debug(workflow)
    assert res.retcode == 0
    for _ in [0, 1]:
        TASK_ID = 1
        res = invoke(
            (
                f"workflow add-task {prj_id} {workflow_id} "
                f"--task-id {TASK_ID} "
                f"--args-non-parallel {args_file}"
            )
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
        f"job submit "
        f"{prj_id} {workflow_id} {dataset_id} "
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
    cmd = f"job show {prj_id} {job_id}"
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
    input_filters_file = str(tmp_path / "input_filters.json")
    args_file = str(tmp_path / "args.json")
    with open(args_file, "w") as f:
        json.dump({"image_dir": "/tmp"}, f)

    with open(input_filters_file, "w") as f:
        json.dump({"raise_error": True}, f)

    res = invoke(
        (
            f"workflow add-task {prj_id} {workflow_id} --task-id {TASK_ID}"
            f" --args-non-parallel {args_file} "
            f"--input-filters {input_filters_file}"
        )
    )
    assert res.retcode == 0
    cmd = f"job submit " f"{prj_id} {workflow_id} {dataset_id}"
    debug(cmd)
    res = invoke(cmd)
    assert res.retcode == 0
    job_id = res.data["id"]

    # Avoid immediately calling `job show` right after `workflow apply`
    time.sleep(1)

    # Verify that status is failed, and that there is a log
    cmd = f"job show {prj_id} {job_id}"
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
    cmd = f"--batch job submit " f"{prj_id} {workflow_id} {dataset_id}"
    debug(cmd)
    res = invoke(cmd)
    assert res.retcode == 0
    debug(res.data)
    assert isinstance(res.data, int)
