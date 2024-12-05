import time
from pathlib import Path
from urllib.request import urlretrieve

import pytest
from devtools import debug

TIMEOUT = 15.0


def test_job_submit(
    invoke,
    project_factory,
    dataset_factory,
    tmp_path: Path,
    testdata_path: Path,
    new_name,
):
    # Collect tasks
    PACKAGE_URL = (
        "https://github.com/fractal-analytics-platform/fractal-server/"
        "raw/main/tests/v2/fractal_tasks_mock/dist/"
        "fractal_tasks_mock-0.0.1-py3-none-any.whl"
    )
    PACKAGE_PATH = "/tmp/fractal_tasks_mock-0.0.1-py3-none-any.whl"
    urlretrieve(PACKAGE_URL, PACKAGE_PATH)

    res = invoke(f"--batch task collect {PACKAGE_PATH} --private")
    assert res.retcode == 0
    activity_id = res.data

    # Create a project
    project = project_factory(name=new_name())
    project_id = project["id"]
    zarr_dir = (tmp_path / "zarr_dir").as_posix()
    dataset = dataset_factory(
        name=new_name(), project_id=project_id, zarr_dir=zarr_dir
    )
    dataset_id = dataset["id"]

    # Wait for task collection to end
    starting_time = time.perf_counter()
    while True:
        res1 = invoke(f"task check-collection {activity_id}")
        if res1.data["status"] == "OK":
            debug(res1.data)
            break
        time.sleep(0.1)
        assert time.perf_counter() - starting_time < TIMEOUT

    wf_json = (testdata_path / "import-export/wf3.json").as_posix()
    res = invoke(
        f"workflow import --project-id {project_id} --json-file {wf_json}"
    )
    workflow = res.data
    workflow_id = workflow["id"]
    debug(workflow)

    FIRST_TASK_INDEX = 0
    LAST_TASK_INDEX = 0
    WORKER_INIT = "export MYVARIABLE=MYVALUE"
    res = invoke(
        f"job submit {project_id} {workflow_id} {dataset_id} "
        f"--start {FIRST_TASK_INDEX} --end {LAST_TASK_INDEX} "
        f'--worker-init "{WORKER_INIT}"'
    )
    assert res.retcode == 0
    job1 = res.data
    job1_id = job1["id"]
    assert job1["status"] == "submitted"
    assert job1["first_task_index"] == FIRST_TASK_INDEX
    assert job1["last_task_index"] == LAST_TASK_INDEX
    assert job1["worker_init"] == WORKER_INIT

    # Check that job completed successfully
    cmd = f"job show {project_id} {job1_id}"
    starting_time = time.perf_counter()
    debug(cmd)
    while True:
        res = invoke(cmd)
        job1 = res.data
        debug(job1)
        assert res.retcode == 0
        if job1["status"] == "done":
            break
        elif job1["status"] == "failed":
            raise RuntimeError(job1)
        time.sleep(0.1)
        assert time.perf_counter() - starting_time < TIMEOUT
    assert job1["log"] is not None

    # Prepare and run a workflow with a failing task
    FIRST_TASK_INDEX = 0
    LAST_TASK_INDEX = 1
    res = invoke(
        f"--batch "
        f"job submit {project_id} {workflow_id} {dataset_id} "
        f"--start {FIRST_TASK_INDEX} --end {LAST_TASK_INDEX} "
        f'--worker-init "{WORKER_INIT}"'
    )
    assert res.retcode == 0
    job2_id = res.data

    # Verify that status is failed, and that there is a log
    cmd = f"--batch  job show {project_id} {job2_id}"
    starting_time = time.perf_counter()
    while True:
        res = invoke(cmd)
        status = res.data
        debug(status)
        assert res.retcode == 0
        if status == "failed":
            break
        time.sleep(0.1)
        assert time.perf_counter() - starting_time < TIMEOUT

    # Run job list with/without --batch
    res = invoke(f"--batch job list {project_id}")
    assert res.retcode == 0
    assert res.data == f"{job1_id} {job2_id}"
    res = invoke(f"job list {project_id}")
    assert res.retcode == 0
    assert set(job["id"] for job in res.data) == set([job1_id, job2_id])

    # Download logs / success
    log1_dir = tmp_path / "log1"
    cmd = (
        f"job download-logs {project_id} {job1_id} "
        f"--output {log1_dir.as_posix()}"
    )
    res = invoke(cmd)
    assert res.retcode == 0
    files = log1_dir.glob("*")
    assert "workflow.log" in [f.name for f in files]

    # Download logs / fail because folder already exists
    log1_dir = tmp_path / "log1"
    cmd = (
        f"job download-logs {project_id} {job1_id} "
        f"--output {log1_dir.as_posix()}"
    )
    res = invoke(cmd)
    assert res.retcode == 1

    # Download logs / fail because of invalid job_id
    cmd = f"job download-logs {project_id} 9999 --output /tmp/invalid/"
    with pytest.raises(SystemExit):
        invoke(cmd)
