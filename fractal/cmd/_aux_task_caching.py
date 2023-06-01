import json
from pathlib import Path
from typing import Any
from typing import Optional

from ..authclient import AuthClient
from ..config import settings
from ..response import check_response

TASKS_CACHE_FILENAME = "tasks"


class FractalCacheError(RuntimeError):
    pass


# Define a useful type
_TaskList = list[dict[str, Any]]


async def _fetch_task_list(client: AuthClient) -> _TaskList:
    """
    Fetch task list through an API request.
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


def _write_task_list(task_list: _TaskList) -> None:
    """
    Write task list to cache file
    """
    cache_dir = Path(f"{settings.FRACTAL_CACHE_PATH}").expanduser()
    cache_dir.mkdir(parents=True, exist_ok=True)
    with (cache_dir / TASKS_CACHE_FILENAME).open("w") as f:
        json.dump(task_list, f, indent=4)


async def refresh_task_cache(client: AuthClient) -> list[dict[str, Any]]:
    """
    Return task list after fetching it, sorting it and writing to cache file.
    """
    task_list = await _fetch_task_list(client)
    task_list = _sort_task_list(task_list)
    _write_task_list(task_list)
    return task_list


def _get_matching_tasks(
    task_list: list[dict[str, Any]],
    *,
    name: str,
    version: Optional[str] = None,
):
    """
    Given a task list, extract all the tasks matching some conditions.
    """

    def _condition(_task):
        if _task["name"] == name:
            if (not version) or (_task["version"] == version):
                return True
            return False

    return [_task for _task in task_list if _condition(_task)]


def _search_in_task_list(
    *,
    task_list: list[dict[str, Any]],
    name: str,
    version: Optional[str] = None,
) -> int:
    """
    Searches for a single task in the `task_list` based on the provided `name`
    and `version`.
    If `version` is not provided, the maximum available version is used.
    Returns the `id` of the single matching task, if found.
    If the task is not found or is not unique, raises a `FractalCacheError`.
    """
    matching_task_list = _get_matching_tasks(
        task_list, name=name, version=version
    )

    if len(matching_task_list) == 0:
        raise FractalCacheError(
            f"There is no task with (name, version)=({name},{version}) "
            "in the cache"
        )
    elif len(matching_task_list) == 1:
        return matching_task_list[0]["id"]
    else:  # i.e. len(matching_task_list) > 1
        if version is None:
            if any(task["version"] is None for task in matching_task_list):
                raise FractalCacheError(
                    "Cannot determine the maximum version in this list "
                    f"(there are one or more None)\n{matching_task_list}"
                )
            max_version = max(matching_task_list, key=lambda x: x["version"])[
                "version"
            ]
            max_version_tasks = []
            for task in task_list:
                if task["version"] == max_version:
                    max_version_tasks.append(task)
            if len(max_version_tasks) == 1:
                return max_version_tasks[0]["id"]
            else:
                raise FractalCacheError(
                    f"Multiple tasks with maximum version {max_version}\n"
                    f"{max_version_tasks}"
                )
        else:
            raise FractalCacheError(
                f"Multiple tasks with version {version}\n{task_list}"
            )


async def get_task_id_from_cache(
    client: AuthClient, task_name: str, version: Optional[str] = None
) -> int:
    """
    Retrieves the `id` of a task from the cache based on the provided
    `task_name` and `version`.
    If `version` is not provided, the maximum available version is used.
    Returns the `id` of the single matching task, if found.
    If the task is not found or is not unique, raises a `FractalCacheError`.
    """

    # Set paths
    cache_dir = Path(f"{settings.FRACTAL_CACHE_PATH}").expanduser()
    cache_file = cache_dir / TASKS_CACHE_FILENAME
    # If cache is missing, create it
    if not cache_file.exists():
        task_list = await refresh_task_cache(client)
    else:
        with cache_file.open("r") as f:
            task_list = json.load(f)
    try:
        task_id = _search_in_task_list(
            task_list=task_list,
            name=task_name,
            version=version,
        )
    except FractalCacheError:
        task_list = await refresh_task_cache(client)
        task_id = _search_in_task_list(
            task_list=task_list,
            name=task_name,
            version=version,
        )
    return task_id
