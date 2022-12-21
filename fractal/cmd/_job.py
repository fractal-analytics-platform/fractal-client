from ..authclient import AuthClient
from ..common.schemas import ApplyWorkflowRead
from ..config import settings
from ..interface import BaseInterface
from ..interface import PrintInterface
from ..interface import RichJsonInterface
from ..response import check_response


async def job_status(
    client: AuthClient,
    job_id: int,
    batch: bool = False,
    do_not_separate_logs: bool = False,
    **kwargs,
) -> BaseInterface:
    """
    Query the status of a workflow-execution job
    """

    res = await client.get(f"{settings.BASE_URL}/job/{job_id}")
    job = check_response(
        res, expected_status_code=200, coerce=ApplyWorkflowRead
    )
    if batch:
        return PrintInterface(retcode=0, data=job.status)
    else:
        data = job.sanitised_dict()
        if do_not_separate_logs or (job.log is None):
            return RichJsonInterface(retcode=0, data=data)
        else:
            log = data.pop("log")
            extra_lines = ["\nThis is the job log:\n", log]
            return RichJsonInterface(
                retcode=0, data=data, extra_lines=extra_lines
            )


async def job_list(
    client: AuthClient,
    project_id: int,
    batch: bool = False,
    **kwargs,
) -> BaseInterface:

    raise NotImplementedError("job_list not implemented")


async def job_download_logs(
    client: AuthClient,
    project_id: int,
    batch: bool = False,
    **kwargs,
) -> BaseInterface:

    raise NotImplementedError("job_download_logs not implemented")
