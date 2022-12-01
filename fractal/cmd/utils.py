import json
import logging
from pathlib import Path
from typing import List

from ..authclient import AuthClient
from ..config import settings
from ..response import check_response


async def get_cached_task_by_name(name: str, client: AuthClient) -> int:

    # Set paths
    cache_dir = Path(f"{settings.FRACTAL_CACHE_PATH}").expanduser()
    cache_file = cache_dir / "tasks"

    # If cache is missing, create it
    cache_up_to_date = False
    if not cache_file.exists():
        await refresh_task_cache(client)
        cache_up_to_date = True

    # Read cache
    with cache_file.open("r") as f:
        task_cache = json.load(f)

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


async def refresh_task_cache(client: AuthClient, **kwargs) -> List[dict]:

    # Get task_list
    res = await client.get(f"{settings.BASE_URL}/task/")
    task_list = check_response(res, expected_status_code=200)

    # Check that there are no name clashes
    names = [task["name"] for task in task_list]
    if len(set(names)) < len(names):
        msg = (
            "Cannot write task-list cache if task names are not unique "
            "(version-based disambiguation will be added in the future).\n"
            f"Current task list includes: {names}"
        )
        logging.warning(msg)
        return task_list

    # Refresh cache of (name,id) pairs
    task_cache = {}
    for task in task_list:
        task_cache[task["name"]] = task["id"]

    # Set paths and create cache folder (if needed)
    cache_dir = Path(f"{settings.FRACTAL_CACHE_PATH}").expanduser()
    cache_file = cache_dir / "tasks"
    cache_dir.mkdir(parents=True, exist_ok=True)

    # Write task cache
    with cache_file.open("w") as f:
        json.dump(task_cache, f, indent=4)

    return task_list
