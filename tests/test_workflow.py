import json

import pytest  # noqa F401
from devtools import debug


async def test_workflow_new(clear_db, testserver, register_user, invoke):
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


async def test_workflow_list(clear_db, testserver, register_user, invoke):
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


async def test_add_task(
    clear_db,
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
