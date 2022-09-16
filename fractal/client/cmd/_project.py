import logging
from typing import Optional

from ..authclient import AuthClient
from ..config import settings
from ..response import check_response
from fractal.common.models import ProjectCreate


async def project_create(
    client: AuthClient,
    name: str,
    path: str,
    dataset: Optional[str] = None,
    **kwargs,
):
    project = ProjectCreate(
        name=name, project_dir=path, default_dataset_name=dataset
    )
    logging.info(project)
    res = await client.post(
        f"{settings.BASE_URL}/project/",
        json=project.dict(),
    )
    data = check_response(res, expected_status_code=201)
    return data


async def project_list():
    pass


async def porject_new():
    pass
