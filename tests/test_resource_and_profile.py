import json
from typing import Any


def get_resource_data(name: str) -> dict[str, Any]:
    return {
        "type": "local",
        "name": name,
        "jobs_local_dir": "/tmp/jobs_local_dir",
        "tasks_local_dir": "/tmp/tasks_local_dir",
        "jobs_poll_interval": 5,
        "jobs_runner_config": {},
        "tasks_python_config": {
            "default_version": "3.12",
            "versions": {"3.12": "/usr/bin/python3.12"},
        },
        "tasks_pixi_config": {},
    }


def test_new_resource_and_profile(tmp_path, invoke_as_superuser, new_name):
    resource_data_1 = get_resource_data(new_name())
    resource_data_2 = get_resource_data(new_name())
    resource_path_1 = tmp_path / "res1.json"
    resource_path_2 = tmp_path / "res2.json"
    with resource_path_1.open("w") as f:
        json.dump(resource_data_1, f)
    with resource_path_2.open("w") as f:
        json.dump(resource_data_2, f)

    res = invoke_as_superuser(f"resource new {resource_path_1.as_posix()}")
    assert res.retcode == 0
    resource_id_1 = res.data["id"]

    res = invoke_as_superuser(
        f"--batch resource new {resource_path_2.as_posix()}"
    )
    assert res.retcode == 0
    resource_id_2 = int(res.data)

    profile_data_1 = dict(resource_type="local", name=f"profile {new_name()}")
    profile_data_2 = dict(resource_type="local", name=f"profile {new_name()}")
    profile_path_1 = tmp_path / "prof1.json"
    profile_path_2 = tmp_path / "prof2.json"
    with profile_path_1.open("w") as f:
        json.dump(profile_data_1, f)
    with profile_path_2.open("w") as f:
        json.dump(profile_data_2, f)

    res = invoke_as_superuser(
        f"profile new {resource_id_1} {profile_path_1.as_posix()}"
    )
    assert res.retcode == 0
    res = invoke_as_superuser(
        f"--batch profile new {resource_id_2} {profile_path_2.as_posix()}"
    )
    int(res.data)
    assert res.retcode == 0
