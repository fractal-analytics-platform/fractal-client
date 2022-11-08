import json
from pathlib import Path
from typing import List
from typing import Optional

from ..authclient import AuthClient
from ..config import settings
from ..interface import BaseInterface
from ..interface import PrintInterface
from ..interface import RichJsonInterface
from ..response import check_response
from ..schemas import TaskCollectPip
from ..schemas import TaskRead
from ..schemas import TaskUpdate


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
        raise KeyError(f"Task {name} not in {cache_file}\n")
    # Case 3: name is missing but cache may be out of date
    else:
        await refresh_task_cache(client)
        with cache_file.open("r") as f:
            task_cache = json.load(f)
        try:
            return task_cache[name]
        except KeyError as e:
            raise KeyError(f"Task {name} not in {cache_file}\n", str(e))


async def refresh_task_cache(client: AuthClient, **kwargs) -> List[dict]:

    # Get task_list
    res = await client.get(f"{settings.BASE_URL}/task/")
    task_list = check_response(res, expected_status_code=200)

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


async def task_list(client: AuthClient, **kwargs) -> RichJsonInterface:
    task_list = await refresh_task_cache(client=client, **kwargs)
    return RichJsonInterface(retcode=0, data=task_list)


async def task_collect_pip(
    client: AuthClient,
    *,
    package: str,
    version: Optional[str] = None,
    python_version: Optional[str],
    package_extras: Optional[str] = None,
    private: Optional[bool] = False,
    **kwargs,
) -> BaseInterface:
    task_collect = TaskCollectPip(
        package=package,
        version=version,
        python_version=python_version or "3.8",
        package_extras=package_extras,
    )

    res = await client.post(
        f"{settings.BASE_URL}/task/collect/pip/?public={not private}",
        json=task_collect.dict(),
    )
    # TODO check response
    if res.status_code == 201:
        return RichJsonInterface(retcode=0, data=res.json())
    else:
        raise RuntimeError(res.json())


async def task_collection_check(
    client: AuthClient, *, installation_path: Path, **kwargs
) -> BaseInterface:
    res = await client.get(
        f"{settings.BASE_URL}/task/collect/{installation_path}"
    )
    if res.status_code == 200:
        return RichJsonInterface(retcode=0, data=res.json())
    else:
        raise RuntimeError(res.json())


async def task_edit(
    client: AuthClient,
    *,
    task_id_or_name: str,
    **task_update_dict,
) -> BaseInterface:
    task_update = TaskUpdate(**task_update_dict)
    payload = task_update.dict(exclude_unset=True)
    if not payload:
        return PrintInterface(retcode=1, data="Nothing to update")

    try:
        task_id = int(task_id_or_name)
    except ValueError:
        task_id = await get_cached_task_by_name(
            name=task_id_or_name, client=client
        )

    res = await client.patch(
        f"{settings.BASE_URL}/task/{task_id}", json=payload
    )
    new_task = check_response(res, expected_status_code=200, coerce=TaskRead)
    return RichJsonInterface(retcode=0, data=new_task.dict())
