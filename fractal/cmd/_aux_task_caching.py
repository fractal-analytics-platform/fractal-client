import json
from pathlib import Path
from typing import Any
from typing import Optional

from ..authclient import AuthClient
from ..config import settings
from ..response import check_response

TASKS_CACHE_FILENAME = "tasks"


# Define a useful type
_TaskList = list[dict[str, Any]]


def _fetch_task_list(client: AuthClient) -> _TaskList:
    """
    Fetch task list through an API request.
    """
    res = client.get(f"{settings.BASE_URL}/task/")
    task_list = check_response(res, expected_status_code=200)
    return task_list


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


def refresh_task_cache(client: AuthClient) -> list[dict[str, Any]]:
    """
    Return task list after fetching it, sorting it and writing to cache file.
    """
    task_list = _fetch_task_list(client)
    task_list = _sort_task_list(task_list)
    _write_task_list(task_list)
    return task_list


def _get_task_id(
    task_list: list[dict],
    *,
    name: Optional[str] = None,
    _id: Optional[int] = None,
    version: Optional[str] = "___VERSION_DEFAULT___",
    owner: Optional[str] = "___OWNER_DEFAULT___",
):
    """
    Given a task list, extract the ID of the single task matching some
    conditions.

    Note: this function also lets the user find tasks with `version` or `owner`
    that are `None`; this won't be necessarily used in other public
    functions.
    """

    def _condition(_task):
        if name and _task["name"] != name:
            return False
        if version != "___VERSION_DEFAULT___" and _task["version"] != version:
            return False
        if owner != "___OWNER_DEFAULT___" and _task["owner"] != owner:
            return False
        if _id and _task["id"] != _id:
            return False
        return True

    tasks = [_task for _task in task_list if _condition(_task)]
    if len(tasks) == 0:
        raise ValueError("No task matches required attributes")
    elif len(tasks) > 1:
        raise ValueError(
            f"Multiple tasks ({len(tasks)}) match required attributes"
        )
    return tasks[0]["id"]


def get_cached_task_by_name(name: str, client: AuthClient) -> int:

    """ """

    # Set paths
    cache_dir = Path(f"{settings.FRACTAL_CACHE_PATH}").expanduser()
    cache_file = cache_dir / TASKS_CACHE_FILENAME

    # If cache is missing, create it
    cache_up_to_date = False
    if not cache_file.exists():
        refresh_task_cache(client)
        cache_up_to_date = True

    # Read task list and create cache of (name, id) pairs
    # FIMXE: this step will be modified in view of
    # https://github.com/fractal-analytics-platform/fractal/issues/345
    with cache_file.open("r") as f:
        task_list = json.load(f)
    task_cache = {}
    for task in task_list:
        if task["name"] in task_cache.keys():
            raise ValueError("Cannot parse task_list")
        task_cache[task["name"]] = task["id"]

    # Look for name in cache
    # Case 1: name exists in cache
    if name in task_cache.keys():
        return task_cache[name]
    # Case 2: name is missing, and cache was just updated
    elif cache_up_to_date:
        raise KeyError(f'Task "{name}" not in {cache_file}\n')
    # Case 3: name is missing but cache may be out of date
    else:
        refresh_task_cache(client)
        with cache_file.open("r") as f:
            task_cache = json.load(f)
        try:
            return task_cache[name]
        except KeyError as e:
            raise KeyError(f'Task "{name}" not in {cache_file}\n', str(e))
