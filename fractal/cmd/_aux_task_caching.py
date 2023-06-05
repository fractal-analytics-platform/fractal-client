import json
from pathlib import Path
from typing import Any
from typing import Optional

from ..authclient import AuthClient
from ..config import settings
from ..response import check_response

TASKS_CACHE_FILENAME = "tasks"


class FractalCacheError(RuntimeError):
    """
    Custom error raised by functions of this module
    """

    pass


# Define a useful type
_TaskList = list[dict[str, Any]]


async def _fetch_task_list(client: AuthClient) -> _TaskList:
    """
    Fetch task list through an API request, and only select a few relevant
    attributes of each task.
    """
    res = await client.get(f"{settings.BASE_URL}/task/")
    task_list = check_response(res, expected_status_code=200)
    return [
        dict(
            id=task["id"],
            name=task["name"],
            version=task["version"],
            owner=task["owner"],
            source=task["source"],
        )
        for task in task_list
    ]


def _sort_task_list(task_list: _TaskList) -> _TaskList:
    """
    Sort tasks according to their (owner, name, version) attributes.
    """
    new_task_list = sorted(
        task_list,
        key=lambda task: (
            task["owner"] or "",
            task["name"],
            task["version"] or "",
        ),
    )
    return new_task_list


def _write_task_list_to_cache(task_list: _TaskList) -> None:
    """
    Write task list to cache file
    """
    cache_dir = Path(f"{settings.FRACTAL_CACHE_PATH}").expanduser()
    cache_dir.mkdir(parents=True, exist_ok=True)
    with (cache_dir / TASKS_CACHE_FILENAME).open("w") as f:
        json.dump(task_list, f, indent=4)


async def refresh_task_cache(client: AuthClient) -> _TaskList:
    """
    Return task list after fetching it, sorting it and writing to cache file.
    """
    task_list = await _fetch_task_list(client)
    task_list = _sort_task_list(task_list)
    _write_task_list_to_cache(task_list)
    return task_list


def _get_matching_tasks(
    task_list: _TaskList,
    *,
    name: str,
    version: Optional[str] = None,
) -> _TaskList:
    """
    Given a task list, extract all the tasks matching some conditions.
    """

    def _condition(_task):
        if _task["name"] == name:
            if (version is None) or (_task["version"] == version):
                return True
            return False
        else:
            return False

    return [_task for _task in task_list if _condition(_task)]


def _format_task_list(task_list: _TaskList) -> str:
    """
    Helper function to print a formatted task list
    """
    header = "  ID, Name, Version, Owner, Source\n"
    formatted_task_list = (
        header
        + "\n".join(
            [
                f'  {task["id"]}, "{task["name"]}", {task["version"]}, {task.get("owner")}, {task["source"]}'  # noqa
                for task in task_list
            ]
        )
        + "\n"
    )
    return formatted_task_list


def _search_in_task_list(
    *,
    task_list: _TaskList,
    name: str,
    version: Optional[str] = None,
) -> int:
    """
    Search for a single task in `task_list` based on the provided `name`
    and `version`, and return its `id`.

    If `version` is not provided, use the maximum available version (that is,
    the latest version).

    If the task is not found or is not unique, raise a `FractalCacheError`.
    """
    matching_task_list = _get_matching_tasks(
        task_list, name=name, version=version
    )
    formatted_matching_task_list = _format_task_list(matching_task_list)

    if len(matching_task_list) == 0:
        formatted_task_list = _format_task_list(task_list)
        if version is not None:
            raise FractalCacheError(
                f'There is no task with (name, version)=("{name}",{version}) '
                f"in the following task list:\n{formatted_task_list}\n"
            )
        else:
            raise FractalCacheError(
                f'There is no task with name "{name}" '
                f"in the following task list:\n{formatted_task_list}\n"
            )
    elif len(matching_task_list) == 1:
        return matching_task_list[0]["id"]
    else:  # i.e. len(matching_task_list) > 1
        if version is not None:
            raise FractalCacheError(
                f"Multiple tasks with version {version} in the following "
                f"task list:\n{formatted_matching_task_list}"
                "Please make your request more specific.\n"
            )
        else:  # i.e. version is None
            if any(task["version"] is None for task in matching_task_list):
                raise FractalCacheError(
                    "Cannot determine the latest version in the following "
                    f"task list:\n{formatted_matching_task_list}"
                    "Please make your request more specific.\n"
                )
            max_version = max(_task["version"] for _task in matching_task_list)
            max_version_tasks = [
                _task
                for _task in matching_task_list
                if _task["version"] == max_version
            ]
            formatted_matching_task_list = _format_task_list(max_version_tasks)
            if len(max_version_tasks) == 1:
                return max_version_tasks[0]["id"]
            else:
                raise FractalCacheError(
                    "Multiple tasks with latest version "
                    f"({max_version}) in the following task "
                    f"list:\n{formatted_matching_task_list}"
                    "Please make your request more specific.\n"
                )


async def get_task_id_from_cache(
    client: AuthClient, task_name: str, version: Optional[str] = None
) -> int:
    """
    Retrieve the `id` of a task from the cache based on the provided
    `task_name` and `version`.

    If `version` is not provided, the latest (i.e. maximum) available version
    is used.

    Return the `id` of the single matching task, if found.

    If the task is not found or is not unique, re-try after refreshing the
    cache, and then raise a `FractalCacheError`.
    """

    # If cache is missing, create it
    cache_dir = Path(f"{settings.FRACTAL_CACHE_PATH}").expanduser()
    cache_file = cache_dir / TASKS_CACHE_FILENAME
    if cache_file.exists():
        with cache_file.open("r") as f:
            task_list = json.load(f)
        already_refreshed_cache = False
    else:
        task_list = await refresh_task_cache(client)
        already_refreshed_cache = True

    try:
        task_id = _search_in_task_list(
            task_list=task_list,
            name=task_name,
            version=version,
        )
    except FractalCacheError as e:
        if already_refreshed_cache:
            # Cache is already up to date, fail
            raise e
        else:
            # Cache may be out-of-date, refresh it and try again
            task_list = await refresh_task_cache(client)
            task_id = _search_in_task_list(
                task_list=task_list,
                name=task_name,
                version=version,
            )
    return task_id
