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
    task_list: list[dict], *, name: str, version: Optional[str] = None
):
    """
    Given a task list, extract all the tasks matching some conditions.
    """

    def _condition(_task):
        if _task["name"] != name:
            if (not version) or (_task["version"] == version):
                return True
            return False

    return [_task for _task in task_list if _condition(_task)]


async def _search_in_task_list(
    *,
    client: AuthClient,
    task_list,
    name: str,
    version: Optional[str] = None,
    allow_cache_refresh: bool = True,
) -> int:

    task_list = _get_matching_tasks(task_list, name=name, version=version)

    if len(task_list) == 0:
        if allow_cache_refresh:
            task_list = await refresh_task_cache(client)
            return await _search_in_task_list(
                client=client,
                task_list=task_list,
                name=name,
                version=version,
                allow_cache_refresh=False,
            )
        else:
            raise ValueError(
                f"There is no task with (name, version)=({name},{version}) "
                "in the cache"
            )
    elif len(task_list) == 1:
        return task_list[0]["id"]
    else:
        if version is None:
            for task in task_list:
                if task["version"] is None:
                    raise ValueError(
                        "Cannot determine max version among this list\n"
                        f"{task_list}"
                    )
            max_version = max(task_list, key=lambda x: x["version"])["version"]
            return await _search_in_task_list(
                client=client,
                task_list=task_list,
                name=name,
                version=max_version,
                allow_cache_refresh=False,
            )
        else:
            raise ValueError(
                "Cannot find a unique task with "
                f"name={name}, version={version} among this list\n"
                f"{task_list}"
            )


async def get_task_id_from_cache(
    client: AuthClient, task_id_or_name: str, version: Optional[str] = None
):
    if str(task_id_or_name).isdigit():
        _id = int(task_id_or_name)
        if version:
            raise Exception("---")  # FIXME
        return _id
    else:
        # Set paths
        cache_dir = Path(f"{settings.FRACTAL_CACHE_PATH}").expanduser()
        cache_file = cache_dir / TASKS_CACHE_FILENAME
        # If cache is missing, create it
        if not cache_file.exists():
            task_list = await refresh_task_cache(client)

        task_id = _search_in_task_list(
            client=client,
            task_list=task_list,
            name=task_id_or_name,
            version=version,
            allow_cache_refresh=True,
        )
        return task_id
