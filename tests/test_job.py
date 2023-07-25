from pathlib import Path

import pytest  # noqa F401
from devtools import debug


LOG = "Here are some logs"


@pytest.mark.parametrize("status", ["done", "failed"])
async def test_job_show(
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
    res = await invoke("project new prj0")
    project_id = res.data["id"]
    wf = await workflow_factory(project_id=project_id)
    debug(wf)
    log = LOG
    workflow_path = tmp_path / f"workflow_{wf.id}"
    workflow_path.mkdir()
    job = await job_factory(
        working_dir=workflow_path.as_posix(),
        worfklow_id=wf.id,
        status=status,
        log=log,
        workflow_dump={
            "name": "my workflow",
            "id": 1,
            "project_id": 1,
            "task_list": [],
        },
    )
    debug(job)

    # Check `job show` output
    cmd = f"job show {project_id} {job.id}"
    debug(cmd)
    res = await invoke(cmd)
    assert res.retcode == 0
    assert res.data["status"] == status
    debug(res.data)
    assert "log" not in res.data
    assert LOG in res.extra_lines
    res.show()

    # Check `job show` output with --batch
    cmd = f"--batch job show {project_id} {job.id}"
    res = await invoke(cmd)
    assert res.retcode == 0
    assert res.data == status

    # Check `job show` output with `--do-not-separate-logs`
    cmd = f"job show {project_id} {job.id} --do-not-separate-logs"
    res = await invoke(cmd)
    debug(res.data)
    assert res.retcode == 0
    assert res.data["status"] == status
    debug(res.data)
    assert res.data["log"] == LOG


async def test_job_list(
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
    res = await invoke("project new prj0")
    project_id = res.data["id"]
    wf_1 = await workflow_factory(project_id=project_id)
    wf_2 = await workflow_factory(project_id=project_id)
    wd_1 = tmp_path / f"workflow_{wf_1.id}"
    wd_2 = tmp_path / f"workflow_{wf_2.id}"
    job1 = await job_factory(
        working_dir=wd_1.as_posix(),
        worfklow_id=wf_1.id,
        status="running",
        workflow_dump={
            "name": "my workflow",
            "id": 1,
            "project_id": 1,
            "task_list": [],
        },
    )
    job2 = await job_factory(
        working_dir=wd_2.as_posix(),
        worfklow_id=wf_2.id,
        status="done",
        workflow_dump={
            "name": "my workflow",
            "id": 1,
            "project_id": 1,
            "task_list": [],
        },
    )
    debug(job1)
    debug(job2)

    # Check `job list` output with --batch option
    cmd = f"--batch job list {project_id}"
    debug(cmd)
    res = await invoke(cmd)
    assert res.retcode == 0
    debug(res.data)
    job_ids = [int(i) for i in res.data.split()]
    assert job1.id in job_ids
    assert job2.id in job_ids

    # Check `job list` output
    cmd = f"job list {project_id}"
    debug(cmd)
    res = await invoke(cmd)
    assert res.retcode == 0
    debug(res.data)
    # There is not much to assert here, apart from successful invocation of the
    # command. We add a res.show() for when pytest is run with the -s flag
    res.show()


async def test_job_download_logs(
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
    res = await invoke("project new prj0")
    project_id = res.data["id"]
    wf = await workflow_factory(project_id=project_id)
    wd = tmp_path / f"workflow_{wf.id}"
    job = await job_factory(
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
    res = await invoke(cmd)
    assert res.retcode == 1

    # Check standard `job download-logs` output
    output = tmp_path / "output_dir_for_logs"
    cmd = (
        f"job download-logs {job.project_id} {job.id} "
        f"--output {str(output)}"
    )
    debug(cmd)
    res = await invoke(cmd)
    assert res.retcode == 0
    debug(res.data)

    # Check that the logfile is there and has the right content
    logfile = output / LOGFILE
    assert logfile.exists()
    with logfile.open("r") as f:
        contents = f.read()
    assert contents == LOG


async def test_job_stop(
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
    res = await invoke("project new prj0")
    project_id = res.data["id"]
    wf = await workflow_factory(project_id=project_id)
    workflow_path = tmp_path / f"workflow_{wf.id}"
    workflow_path.mkdir()
    job = await job_factory(
        working_dir=workflow_path.as_posix(),
        worfklow_id=wf.id,
    )
    debug(job)

    # Call `job stop` (this will fail because FRACTAL_RUNNER_BACKEND="local"
    cmd = f"job stop {project_id} {job.id}"
    debug(cmd)
    with pytest.raises(SystemExit):
        res = await invoke(cmd)
    debug(caplog.text)
    assert "Stopping a job execution is not implemented" in caplog.text
    caplog.clear()
