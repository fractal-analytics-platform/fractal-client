from typing import Any


def _trim_dict(raw_data, allowed_keys):
    if raw_data:
        return {key: raw_data[key] for key in allowed_keys}
    else:
        return raw_data


def _simplify_project(raw_data: dict[str, Any]) -> dict[str, Any]:
    return _trim_dict(
        raw_data=raw_data, allowed_keys=["name", "id", "timestamp_created"]
    )


def _simplify_resource(raw_data: dict[str, Any]) -> dict[str, Any]:
    return _trim_dict(
        raw_data=raw_data, allowed_keys=["id", "dataset_id", "path"]
    )


def _simplify_dataset(raw_data: dict[str, Any]) -> dict[str, Any]:
    from devtools import debug

    debug("BBBBB", raw_data)
    tmp = _trim_dict(
        raw_data=raw_data,
        allowed_keys=[
            "id",
            "name",
            "type",
            "read_only",
            "project_id",
            "timestamp_created",
        ],
    )
    tmp["resource_list"] = _simplify_resource(raw_data["resource_list"])
    return tmp


def _simplify_job(raw_data: dict[str, Any]) -> dict[str, Any]:
    return _trim_dict(
        raw_data=raw_data,
        allowed_keys=[
            "id",
            "workflow_id",
            "project_id",
            "input_dataset_id",
            "status",
            "start_timestamp",
            "end_timestamp",
        ],
    )


def _simplify_task(raw_data: dict[str, Any]) -> dict[str, Any]:
    return _trim_dict(
        raw_data=raw_data,
        allowed_keys=["name", "command", "source", "version", "id"],
    )


def _simplify_wftask(raw_data: dict[str, Any]) -> dict[str, Any]:
    tmp = _trim_dict(raw_data=raw_data, allowed_keys=["task_id", "id"])
    tmp["task"] = _simplify_task(raw_data["task"])
    return tmp


def _simplify_workflow(raw_data: dict[str, Any]) -> dict[str, Any]:
    tmp = _trim_dict(
        raw_data=raw_data,
        allowed_keys=["id", "project_id", "timestamp_created"],
    )
    tmp["task_list"] = [
        _simplify_wftask(wftask) for wftask in raw_data["task_list"]
    ]
    return tmp
