import json
from typing import Any
from typing import Dict
from typing import Optional

from ..authclient import AuthClient
from ..config import settings
from ..interface import BaseInterface
from ..interface import PrintInterface
from ..interface import RichJsonInterface
from ..response import check_response
from fractal.common.models import DatasetRead
from fractal.common.models import DatasetUpdate
from fractal.common.models import ResourceCreate
from fractal.common.models import ResourceRead


async def dataset_show(
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


async def dataset_edit(
    client: AuthClient,
    *,
    project_id: int,
    dataset_id: int,
    dataset_update_dict: Dict[str, Any],
) -> BaseInterface:

    metadata_filename = dataset_update_dict.get("metadata")
    if metadata_filename == "none":
        dataset_update_dict.update(meta={})
    elif metadata_filename is not None:
        meta = json.loads(metadata_filename)
        dataset_update_dict.update(meta=meta)

    dataset_update = DatasetUpdate(**dataset_update_dict)
    payload = dataset_update.dict(exclude_unset=True)
    if not payload:
        return PrintInterface(retcode=1, output="Nothing to update")

    res = await client.patch(
        f"{settings.BASE_URL}/project/{project_id}/{dataset_id}",
        json=dataset_update.dict(exclude_unset=True),
    )
    new_dataset = check_response(
        res, expected_status_code=200, coerce=DatasetRead
    )
    return RichJsonInterface(retcode=0, data=new_dataset.dict())
