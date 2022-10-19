import logging
import os

from sqlalchemy.ext.asyncio import AsyncSession

from ... import __VERSION__
from ...config_runner import settings
from ..models import Dataset
from ..models import Workflow
from .common import auto_output_dataset  # noqa: F401
from .common import close_job_logger
from .common import set_job_logger
from .common import validate_workflow_compatibility  # noqa: F401


if settings.RUNNER_BACKEND == "process":
    from .process import process_workflow
else:

    def no_function(*args, **kwarsg):
        raise NotImplementedError(
            f"Runner backend {settings.RUNNER_BACKEND} not implemented"
        )

    process_workflow = no_function


async def submit_workflow(
    *,
    db: AsyncSession,
    workflow: Workflow,
    input_dataset: Dataset,
    output_dataset: Dataset,
    job_id: int,
    username: str = None,
    worker_init: str = None,
):
    """
    Prepares a workflow and applies it to a dataset

    Arguments
    ---------
    db: (AsyncSession):
        Asynchronous database session
    output_dataset (Dataset | str) :
        the destination dataset of the workflow. If not provided, overwriting
        of the input dataset is implied and an error is raised if the dataset
        is in read only mode. If a string is passed and the dataset does not
        exist, a new dataset with that name is created and within it a new
        resource with the same name.
    """

    input_paths = input_dataset.paths
    output_path = output_dataset.paths[0]

    workflow_id = workflow.id

    WORKFLOW_DIR = settings.RUNNER_ROOT_DIR / "workflow_{workflow_id:06d}"
    if not os.path.isdir(WORKFLOW_DIR):
        os.mkdir(WORKFLOW_DIR)

    logger = set_job_logger(
        logger_name=f"WF{workflow_id}",
        log_file_path=WORKFLOW_DIR / "workflow.log",
        level=logging.INFO,
        formatter=logging.Formatter("%(asctime)s; %(levelname)s; %(message)s"),
    )

    logger.info(f"fractal_server.__VERSION__: {__VERSION__}")
    logger.info(f"START workflow {workflow.name}")
    output_dataset.meta = await process_workflow(
        workflow=workflow,
        input_paths=input_paths,
        output_path=output_path,
        input_metadata=input_dataset.meta,
        username=username,
        logger=logger,
        workflow_dir=WORKFLOW_DIR,
    )

    logger.info(f'END workflow "{workflow.name}"')
    close_job_logger(logger)
    db.add(output_dataset)

    await db.commit()
