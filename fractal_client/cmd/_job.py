import json
import logging
import os
from pathlib import Path
from zipfile import ZipFile

from ..authclient import AuthClient
from ..interface import Interface
from ..response import check_response


def get_job(
    client: AuthClient,
    *,
    project_id: int,
    job_id: int,
    batch: bool = False,
) -> Interface:
    """
    Query the status of a workflow-execution job
    """

    res = client.get(f"api/v2/project/{project_id}/job/{job_id}/")
    job = check_response(res, expected_status_code=200)
    if batch:
        return Interface(retcode=0, data=job["status"])
    else:
        return Interface(retcode=0, data=job)


def get_job_list(
    client: AuthClient, *, project_id: int, batch: bool = False
) -> Interface:

    res = client.get(f"api/v2/project/{project_id}/job/")
    jobs = check_response(res, expected_status_code=200)

    if batch:
        job_ids = " ".join(str(job["id"]) for job in jobs)
        return Interface(retcode=0, data=job_ids)
    else:
        return Interface(retcode=0, data=jobs)


def get_job_logs(
    client: AuthClient,
    *,
    project_id: int,
    job_id: int,
    output_folder: str,
) -> Interface:

    # Check that output_folder does not already exist
    if Path(output_folder).exists():
        return Interface(
            retcode=1, data=f"ERROR: {output_folder=} already exists"
        )

    # Send request to server
    res = client.get(f"api/v2/project/{project_id}/job/{job_id}/download/")

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
    zipped_archive_path = output_folder + "_tmp.zip"
    with open(zipped_archive_path, "wb") as f:
        f.write(res.content)

    # Unzip the log archive
    unzipped_archived_path = output_folder
    os.mkdir(unzipped_archived_path)
    with ZipFile(zipped_archive_path, mode="r") as zipfile:
        zipfile.extractall(path=unzipped_archived_path)

    # Remove zipped temporary file
    os.unlink(zipped_archive_path)

    return Interface(retcode=0, data=f"Logs downloaded to {output_folder=}")


def stop_job(client: AuthClient, *, project_id: int, job_id: int) -> Interface:
    """
    Stop a workflow-execution job
    """

    res = client.get(f"api/v2/project/{project_id}/job/{job_id}/stop/")
    check_response(res, expected_status_code=202)
    return Interface(
        retcode=0, data="Correctly called the job-stopping endpoint"
    )


def job_submit(
    client: AuthClient,
    *,
    project_id: int,
    workflow_id: int,
    dataset_id: int,
    first_task_index: int | None = None,
    last_task_index: int | None = None,
    worker_init: str | None = None,
    attribute_filters_json: str | None = None,
    type_filters_json: str | None = None,
    batch: bool = False,
) -> Interface:

    job_submit = dict()
    # Prepare JobV2 object, without None attributes
    if worker_init is not None:
        job_submit["worker_init"] = worker_init
    if first_task_index is not None:
        job_submit["first_task_index"] = first_task_index
    if last_task_index is not None:
        job_submit["last_task_index"] = last_task_index

    if attribute_filters_json is not None:
        with Path(attribute_filters_json).open("r") as f:
            job_submit["attribute_filters"] = json.load(f)
    if type_filters_json is not None:
        with Path(type_filters_json).open("r") as f:
            job_submit["type_filters"] = json.load(f)

    # Prepare query parameters
    query_parameters = f"workflow_id={workflow_id}" f"&dataset_id={dataset_id}"

    res = client.post(
        (f"api/v2/project/{project_id}/job/" f"submit/?{query_parameters}"),
        json=job_submit,
    )
    job_read = check_response(res, expected_status_code=202)

    if batch:
        return Interface(retcode=0, data=job_read["id"])
    else:
        return Interface(retcode=0, data=job_read)
