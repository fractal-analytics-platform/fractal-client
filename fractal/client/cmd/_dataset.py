from typing import Optional

from ..authclient import AuthClient
from ..config import settings
from ..interface import BaseInterface
from ..interface import RichJsonInterface
from ..response import check_response
from fractal.common.models import ResourceCreate
from fractal.common.models import ResourceRead


async def dataset_show(
    client: AuthClient,
    project_id: int,
    dataset_id: int,
    **kwargs,
) -> BaseInterface:
    raise NotImplementedError


async def dataset_edit(
    client: AuthClient,
    project_id: int,
    dataset_id: int,
    **kwargs,
) -> BaseInterface:
    raise NotImplementedError


async def dataset_add_resource(
    client: AuthClient,
    *,
    project_id: int,
    dataset_id: int,
    path: str,
    glob_pattern: Optional[str] = "",
    **kwargs,
) -> BaseInterface:
    resource = ResourceCreate(path=path, glob_pattern=glob_pattern)

    res = await client.post(
        f"{settings.BASE_URL}/project/{project_id}/{dataset_id}",
        json=resource.dict(),
    )
    new_resource = check_response(
        res, expected_status_code=201, coerce=ResourceRead
    )
    return RichJsonInterface(retcode=0, data=new_resource.dict())
