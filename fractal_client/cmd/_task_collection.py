import json
import logging
import sys
from pathlib import Path

from fractal_client.authclient import AuthClient
from fractal_client.interface import Interface
from fractal_client.response import check_response


def task_collect_pip(
    client: AuthClient,
    *,
    package: str,
    package_version: str | None = None,
    python_version: str | None = None,
    package_extras: str | None = None,
    pinned_dependency: list[str] | None = None,
    private: bool = False,
    batch: bool = False,
) -> Interface:
    # Construct TaskCollectPip object
    task_collect = dict()
    if package_version:
        task_collect["package_version"] = package_version
    if python_version:
        task_collect["python_version"] = python_version
    if package_extras:
        task_collect["package_extras"] = package_extras
    if pinned_dependency:
        for pin in pinned_dependency:
            if len(pin.split("=")) != 2:
                logging.error(
                    f"Invalid pin: {pin}.\nPins must be written as "
                    "'--pinned-dependency PACKAGE_NAME=PACKAGE_VERSION'"
                )
                sys.exit(1)
        task_collect["pinned_package_versions"] = json.dumps(
            {
                _name: _version
                for _name, _version in (
                    p.split("=") for p in pinned_dependency
                )
            }
        )

    is_private = "?private=true" if private else ""
    endpoint_url = f"api/v2/task/collect/pip/{is_private}"
    if package.endswith(".whl"):
        with open(package, "rb") as f:
            file = {
                "file": (
                    Path(package).name,
                    f.read(),
                    "application/zip",
                )
            }
        res = client.post(
            endpoint_url,
            data=task_collect,
            files=file,
        )
    else:
        task_collect["package"] = package
        res = client.post(
            endpoint_url,
            data=task_collect,
        )
    task_group_activity = check_response(res, expected_status_code=202)
    if batch:
        return Interface(retcode=0, data=task_group_activity["id"])
    else:
        return Interface(retcode=0, data=task_group_activity)


def task_collect_custom(
    client: AuthClient,
    *,
    label: str,
    python_interpreter: str,
    manifest: str,
    version: str | None = None,
    package_name: str | None = None,
    package_root: str | None = None,
    private: bool = False,
    batch: bool = False,
) -> Interface:
    try:
        with open(manifest) as f:
            manifest_dict = json.load(f)
    except FileNotFoundError as e:
        raise FileNotFoundError(
            f"Fractal Client cannot find the file {manifest}. "
            "Note that the file must be on the same machine where Fractal "
            f"Client is running.\nOriginal error: {e}."
        )

    task_collect = dict(
        label=label,
        python_interpreter=python_interpreter,
        manifest=manifest_dict,
    )
    if version:
        task_collect["version"] = version
    if package_name:
        task_collect["package_name"] = package_name
    if package_root:
        task_collect["package_root"] = package_root
    is_private = "?private=true" if private else ""

    res = client.post(
        f"api/v2/task/collect/custom/{is_private}",
        json=task_collect,
    )

    task_list = check_response(
        res, expected_status_code=201, redact_long_payload=True
    )

    if batch:
        task_ids = [str(task["id"]) for task in task_list]
        return Interface(retcode=0, data=" ".join(task_ids))
    else:
        return Interface(retcode=0, data=task_list)


def show_task_group_activity(
    client: AuthClient,
    *,
    task_group_activity_id: int,
    include_logs: bool,
) -> Interface:
    res = client.get(f"api/v2/task-group/activity/{task_group_activity_id}/")
    task_group_activity = check_response(res, expected_status_code=200)

    # Remove key-value pairs with None value
    if include_logs is False:
        task_group_activity["log"] = None

    return Interface(retcode=0, data=task_group_activity)
