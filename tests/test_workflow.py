import pytest  # noqa F401


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
