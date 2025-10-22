import json
from typing import Any

from devtools import debug


def get_resource(name: str) -> dict[str, Any]:
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


def test_register_resource(tmp_path, invoke_as_superuser):
    resource = get_resource("res1")
    resource_path = tmp_path / "res1.json"
    with resource_path.open("w") as f:
        json.dump(resource, f)

    res = invoke_as_superuser(f"resource new {resource_path.as_posix()}")
    debug(res.data)
    assert res.retcode == 0
