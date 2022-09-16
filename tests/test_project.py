from devtools import debug


async def test_project_create(
    clear_db, testserver, register_user, clisplit, invoke
):
    PROJECT_NAME = "project_name"
    PROJECT_PATH = "project_path"
    res = await invoke(f"project new {PROJECT_NAME} {PROJECT_PATH}")
    debug(res)
    assert res.data["name"] == PROJECT_NAME
    assert res.data["project_dir"] == PROJECT_PATH


async def test_project_create_batch(
    clear_db, testserver, register_user, clisplit, invoke
):
    res = await invoke("--batch project new project_name project_path")
    debug(res)
    debug(res.output)
    assert res.output == "1"


async def test_project_list(
    clear_db, testserver, register_user, clisplit, invoke
):
    res = await invoke("project list")
    debug(res)
    res.show()

    await invoke("--batch project new prj0 prj_path0")
    await invoke("--batch project new prj1 prj_path1")

    res = await invoke("project list")
    debug(res)
    res.show()
    assert False
