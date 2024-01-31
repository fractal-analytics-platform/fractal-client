from typing import Any


def _trim_dict(raw_data, allowed_keys):
    return {key: raw_data[key] for key in allowed_keys}


def _simplify_project(raw_data: dict[str, Any]) -> dict[str, Any]:
    return _trim_dict(
        raw_data=raw_data, allowed_keys=["name", "id", "timestamp_created"]
    )


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
