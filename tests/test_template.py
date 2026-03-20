import json
from pathlib import Path

import pytest

FRACTAL_DEFAULT_GROUP_NAME = "All"


def test_template_commands(
    invoke,
    invoke_as_superuser,
    new_name,
    workflow_factory,
    task_factory,
    tmp_path: Path,
    capsys,
):
    res = invoke_as_superuser("group list --user-ids")
    default_group = next(
        group
        for group in res.data
        if group["name"] == FRACTAL_DEFAULT_GROUP_NAME
    )
    default_group_id = default_group["id"]

    res = invoke(f"project new {new_name()}")
    project = res.data
    workflow = workflow_factory(name=new_name(), project_id=project["id"])
    task = task_factory(name=new_name(), command_parallel="pwd", version="1")

    invoke(
        "workflow add-task "
        f"{project['id']} {workflow['id']} --task-id {task['id']}"
    )

    # Template new (from workflow_id)
    with pytest.raises(SystemExit):
        res = invoke(f"template new --workflow-id {workflow['id']}")
    TEMPLATE_NAME = new_name()
    res = invoke(
        f"template new --workflow-id {workflow['id']} "
        f"--name {TEMPLATE_NAME} --version 1"
    )
    assert res.retcode == 0
    template1_id = res.data["id"]
    res = invoke(
        f"--batch template new --workflow-id {workflow['id']} "
        f"--name {TEMPLATE_NAME} --version 2 "
        f"--user-group-id {default_group_id}"
    )
    assert res.retcode == 0
    template2_id = int(res.data)

    # Template show
    res = invoke(f"template show {template1_id}")
    assert res.retcode == 0
    assert res.data["user_group_id"] is None
    res = invoke(f"template show {template2_id}")
    assert res.retcode == 0
    assert res.data["user_group_id"] == default_group_id

    # Template delete
    res = invoke(f"template delete {template2_id}")
    assert res.retcode == 0
    with pytest.raises(SystemExit):
        invoke(f"template show {template2_id}")

    # Template export
    template_filename = tmp_path / "template.json"
    res = invoke(f"template export {template1_id} {template_filename}")
    assert res.retcode == 0

    # Template new (from JSON file)
    with pytest.raises(SystemExit):
        # Same (user, name, version) of template1
        invoke(
            f"template new --json-file {template_filename} "
            f"--user-group-id {default_group_id}"
        )
    res = invoke(
        "--batch "
        f"template new --json-file {template_filename} --name {new_name()} "
        f"--user-group-id {default_group_id}"
    )
    assert res.retcode == 0
    template3_id = int(res.data)
    res = invoke(
        "--batch "
        f"template new --json-file {template_filename} --version 2 "
        f"--user-group-id {default_group_id}"
    )
    assert res.retcode == 0

    # Import workflow from template
    with pytest.raises(SystemExit):
        # Same (project_id, name) of workflow
        invoke(f"workflow import-from-template {project['id']} {template3_id}")
    res = invoke(f"--batch workflow list {project['id']}")
    assert res.retcode == 0
    assert len(res.data.split()) == 1
    res = invoke(
        f"workflow import-from-template {project['id']} {template3_id} "
        f"--name {new_name()}"
    )
    assert res.retcode == 0
    res = invoke(f"--batch workflow list {project['id']}")
    assert res.retcode == 0
    assert len(res.data.split()) == 2
    res = invoke(
        "--batch "
        f"workflow import-from-template {project['id']} {template3_id} "
        f"--name {new_name()}"
    )
    assert res.retcode == 0
    assert len(res.data.split()) == 2
    res = invoke(f"--batch workflow list {project['id']}")
    assert res.retcode == 0
    assert len(res.data.split()) == 3

    # flexibility
    with template_filename.open("r") as f:
        template_to_import = json.load(f)
    NEW_VERSION = "1234"
    template_to_import["data"]["task_list"][0]["task"]["version"] = NEW_VERSION
    template_to_import["name"] = new_name()
    template_to_import["data"]["name"] = new_name()
    with template_filename.open("w") as f:
        json.dump(template_to_import, f)
    res = invoke(f"template new --json-file {template_filename}")
    assert res.retcode == 0
    template4_id = res.data["id"]
    with pytest.raises(SystemExit):
        invoke(
            f"workflow import-from-template {project['id']} {template4_id} "
            f"--name {new_name()}"
        )
    assert (
        f"Task '{task['name']}' "
        f"(package '{task['name']}', version '{NEW_VERSION}') not available. "
        f"Available versions: ['{task['version']}']."
    ) in capsys.readouterr().out
