import json
from pathlib import Path
from typing import Any

from ..authclient import AuthClient
from ..config import settings
from ..response import check_response

TASKS_CACHE_FILENAME = "tasks"


# Define a useful type
_TaskList = list[dict[str, Any]]


async def _fetch_task_list(client: AuthClient) -> _TaskList:
    """
    Fetch task list through an API request.
    """
    res = await client.get(f"{settings.BASE_URL}/task/")
    task_list = check_response(res, expected_status_code=200)
    return task_list


def _write_task_list(task_list: _TaskList) -> None:
    """
    Write task list to cache file
    """
    cache_dir = Path(f"{settings.FRACTAL_CACHE_PATH}").expanduser()
    cache_dir.mkdir(parents=True, exist_ok=True)
    with (cache_dir / TASKS_CACHE_FILENAME).open("w") as f:
        json.dump(task_list, f, indent=4)


def _sort_task_list(task_list: _TaskList) -> _TaskList:

    new_task_list = sorted(
        task_list,
        key=lambda task: (
            task["owner"] or "",
            task["name"],
            task["version"] or "",
        ),
    )
    return new_task_list


async def refresh_task_cache(client: AuthClient) -> list[dict[str, Any]]:
    """
    Fetch task list, write cache file, return task list.
    """
    task_list = await _fetch_task_list(client)
    task_list = _sort_task_list(task_list)
    _write_task_list(task_list)
    return task_list


async def get_cached_task_by_name(name: str, client: AuthClient) -> int:

    # Set paths
    cache_dir = Path(f"{settings.FRACTAL_CACHE_PATH}").expanduser()
    cache_file = cache_dir / TASKS_CACHE_FILENAME

    # If cache is missing, create it
    cache_up_to_date = False
    if not cache_file.exists():
        await refresh_task_cache(client)
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
        await refresh_task_cache(client)
        with cache_file.open("r") as f:
            task_cache = json.load(f)
        try:
            return task_cache[name]
        except KeyError as e:
            raise KeyError(f'Task "{name}" not in {cache_file}\n', str(e))
