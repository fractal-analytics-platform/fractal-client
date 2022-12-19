from typing import Optional

from ..authclient import AuthClient
from ..common.schemas import StateRead
from ..common.schemas import TaskCollectPip
from ..common.schemas import TaskRead
from ..common.schemas import TaskUpdate
from ..config import settings
from ..interface import BaseInterface
from ..interface import PrintInterface
from ..interface import RichJsonInterface
from ..response import check_response
from .utils import get_cached_task_by_name
from .utils import refresh_task_cache


async def task_list(client: AuthClient, **kwargs) -> RichJsonInterface:
    task_list = await refresh_task_cache(client=client, **kwargs)
    return RichJsonInterface(retcode=0, data=task_list)


async def task_collect_pip(
    client: AuthClient,
    *,
    package: str,
    package_version: Optional[str] = None,
    python_version: Optional[str],
    package_extras: Optional[str] = None,
    private: Optional[bool] = False,
    batch: bool = False,
    **kwargs,
) -> BaseInterface:
    task_collect = TaskCollectPip(
        package=package,
        version=package_version,
        python_version=python_version,
        package_extras=package_extras,
    )

    res = await client.post(
        f"{settings.BASE_URL}/task/collect/pip/?public={not private}",
        json=task_collect.dict(),
    )

    state = check_response(
        res, expected_status_code=[200, 201], coerce=StateRead
    )
    if batch:
        output = f"{state.id} {state.data['venv_path']}"
        return PrintInterface(retcode=0, data=output)
    else:
        return RichJsonInterface(retcode=0, data=state.sanitised_dict())


async def task_collection_check(
    client: AuthClient, *, state_id: int, verbose: bool, **kwargs
) -> BaseInterface:
    res = await client.get(
        f"{settings.BASE_URL}/task/collect/{state_id}?verbose={verbose}"
    )
    state = check_response(res, expected_status_code=200, coerce=StateRead)
    return RichJsonInterface(retcode=0, data=state.sanitised_dict())


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
