import logging
import os
from pathlib import Path
from zipfile import ZipFile

from rich.table import Table

from ..authclient import AuthClient
from ..common.schemas import ApplyWorkflowRead
from ..config import settings
from ..interface import BaseInterface
from ..interface import PrintInterface
from ..interface import RichConsoleInterface
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

    res = await client.get(f"{settings.BASE_URL}/project/{project_id}/jobs/")
    jobs = check_response(res, expected_status_code=200)

    # Coerce to a list of ApplyWorkflowRead objects
    jobs = [ApplyWorkflowRead(**job) for job in jobs]

    if batch:
        job_ids = " ".join(str(job.id) for job in jobs)
        return PrintInterface(retcode=0, data=job_ids)
    else:
        table = Table(title=f"Job list for project {project_id}")
        kwargs = dict(style="white", justify="center")
        table.add_column("id", **kwargs)
        table.add_column("start_timestamp", **kwargs)
        table.add_column("status", **kwargs)
        table.add_column("workflow_id", **kwargs)
        table.add_column("working_dir", **kwargs)

        for j in jobs:
            table.add_row(
                str(j.id),
                j.sanitised_dict()["start_timestamp"],
                j.status,
                str(j.workflow_id),
                j.working_dir,
            )
        return RichConsoleInterface(retcode=0, data=table)


async def job_download_logs(
    client: AuthClient,
    job_id: int,
    output: str,
    **kwargs,
) -> BaseInterface:

    # Check that output folder does not already exist
    if Path(output).exists():
        return PrintInterface(
            retcode=1, data=f"ERROR: {output} already exists"
        )

    # Send request to server
    res = await client.get(f"{settings.BASE_URL}/job/download/{job_id}")

    # NOTE: We cannot use our default check_response here, because res._content
    # is binary. Therefore we check the status code by hand
    if res.status_code != 200:
        logging.error(f"Server returned {res.status_code}")
        logging.error(
            f"Original request: {res._request.method} {res._request.url}"
        )
        logging.error(
            f"Original payload: {res._request._content.decode('utf-8')}"
        )
        logging.error("Terminating.\n")
        exit(1)

    # Check the content-type entry in the response headers
    content_type = res.headers["content-type"]
    expected_content_type = "application/x-zip-compressed"
    if content_type != expected_content_type:
        logging.error(
            f"Unexpected {content_type=} in headers of server "
            f"response, instead of {expected_content_type=}"
        )
        logging.error(
            f"Original request: {res._request.method} {res._request.url}"
        )
        logging.error(
            f"Original payload: {res._request._content.decode('utf-8')}"
        )
        logging.error("Terminating.\n")
        exit(1)

    # Write response into a temporary zipped file
    zipped_archive_path = output + "_tmp.zip"
    with open(zipped_archive_path, "wb") as f:
        f.write(res.content)

    # Unzip the log archive
    unzipped_archived_path = output
    os.mkdir(unzipped_archived_path)
    with ZipFile(zipped_archive_path, mode="r") as zipfile:
        zipfile.extractall(path=unzipped_archived_path)

    # Remove zipped temporary file
    os.unlink(zipped_archive_path)

    return PrintInterface(
        retcode=0, data=f"Logs downloaded to {output} folder"
    )
