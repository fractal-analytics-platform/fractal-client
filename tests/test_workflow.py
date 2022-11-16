import json

import pytest  # noqa F401
from devtools import debug


async def test_workflow_new(testserver, register_user, invoke):
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


async def test_workflow_list(testserver, register_user, invoke):
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
    testserver, register_user, invoke, tmp_path
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
    testserver,
    invoke,
    register_user,
    task_factory,
    workflow_factory,
    tmp_path,
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

    custom_args = dict(custom="args")
    args_file = tmp_path / "args_file.json"
    with args_file.open("w") as f:
        json.dump(custom_args, f)

    cmd = f"workflow add-task {wf.id} {t.id} --args-file {args_file}"
    debug(cmd)
    res = await invoke(cmd)
    assert res.retcode == 0
    debug(res.data)
    assert res.data["task_list"][0]["args"] == custom_args


async def test_edit_workflow_task(
    testserver,
    invoke,
    register_user,
    task_factory,
    workflow_factory,
    tmp_path,
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
    payload = dict(
        args={"some_arg": "some_value"}, meta={"executor": "cpu-low"}
    )

    json_file = tmp_path / "payload.json"
    with json_file.open("w") as f:
        json.dump(payload, f)

    # Edit workflow task
    debug(res.data)
    workflow_task_id = res.data["task_list"][0]["id"]
    cmd = (
        f"workflow edit-task {wf.id} {workflow_task_id} "
        f"--json-file {json_file}"
    )
    debug(cmd)
    res = await invoke(cmd)
    assert res.retcode == 0
    assert res.data["args"] == payload["args"]
    assert res.data["meta"] == payload["meta"]

    # Check that also the workflow in the db was correctly updated
    res = await invoke(f"workflow show {wf.id}")
    assert res.retcode == 0
    debug(res.data)
    assert res.data["task_list"][0]["args"] == payload["args"]
    assert res.data["task_list"][0]["meta"] == payload["meta"]

    # Check if the correct error is raised where parallelization_level
    # is set
    payload_error = dict(meta={"parallelization_level": "XXX"})

    json_file = tmp_path / "payload_error.json"
    with json_file.open("w") as f:
        json.dump(payload_error, f)

    workflow_task_id = res.data["task_list"][0]["id"]
    cmd = (
        f"workflow edit-task {wf.id} {workflow_task_id} "
        f"--json-file {json_file}"
    )
    debug(cmd)
    with pytest.raises(ValueError):
        res = await invoke(cmd)
