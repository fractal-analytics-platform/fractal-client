import json
import time
from pathlib import Path
from urllib.request import urlretrieve

import pytest
from devtools import debug
from httpx import Request
from httpx import Response

TIMEOUT = 25.0


def test_job_submit(
    invoke,
    project_factory,
    tmp_path: Path,
    testdata_path: Path,
    new_name,
    invoke_as_superuser,
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

    type_filters = {"a": True, "b": False}
    type_filters_file = tmp_path / "type_filters.json"
    with type_filters_file.open("w") as f:
        json.dump(type_filters, f)

    # Add a `tmp_path` to the user's `project_dirs`
    res = invoke("--batch user whoami")
    user_id = int(res.data)
    res = invoke_as_superuser(
        f"user edit {user_id} --add-project-dir {tmp_path.as_posix()}"
    )
    assert res.retcode == 0

    # Create dataset (batch)
    res = invoke(
        "--batch "
        f"project add-dataset {project_id} {new_name()} "
        f"--project-dir {tmp_path.as_posix()} --zarr-subfolder zarr"
    )
    assert res.retcode == 0
    dataset_id = int(res.data)

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
    cmd = f"--batch job show {project_id} {job2_id}"
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
    assert res.data == f"{job2_id} {job1_id}"
    res = invoke(f"job list {project_id}")
    assert res.retcode == 0
    assert {job["id"] for job in res.data} == {job1_id, job2_id}

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

    # --attribute-filters-json and --type-filters-json
    attribute_filters = {"x": [1, 2], "y": ["foo", "bar"]}
    attribute_filters_file = tmp_path / "attribute_filters.json"
    with attribute_filters_file.open("w") as f:
        json.dump(attribute_filters, f)

    type_filters = {"x": True, "y": False}
    type_filters_file = tmp_path / "type_filters.json"
    with type_filters_file.open("w") as f:
        json.dump(type_filters, f)

    res = invoke(
        f"job submit {project_id} {workflow_id} {dataset_id} "
        f"--attribute-filters-json {attribute_filters_file} "
        f"--type-filters-json {type_filters_file}"
    )
    assert res.retcode == 0
    assert res.data["attribute_filters"] == attribute_filters
    assert res.data["type_filters"] == type_filters


def test_job_stop(invoke, caplog):
    with pytest.raises(SystemExit):
        invoke("job stop 123456 1234546")
    EXPECTED_MSG = (
        "Stopping a job execution is not implemented "
        "for FRACTAL_RUNNER_BACKEND=local"
    )
    assert EXPECTED_MSG in caplog.text


def test_job_logs_wrong_content_type(monkeypatch, invoke, tmp_path, caplog):
    def patched_get(*args, **kwargs):
        request = Request(method="GET", url="fake-url")
        response = Response(
            status_code=200,
            headers={"content-type": "text/html"},
            request=request,
        )
        debug(response)
        return response

    import fractal_client.cmd._job

    monkeypatch.setattr(fractal_client.cmd._job.AuthClient, "get", patched_get)

    outdir = (tmp_path / "logs").as_posix()
    caplog.clear()
    with pytest.raises(SystemExit):
        invoke(f"job download-logs 111111111 222222222 --output {outdir}")
    debug(caplog.text)
    assert "Unexpected content_type" in caplog.text
