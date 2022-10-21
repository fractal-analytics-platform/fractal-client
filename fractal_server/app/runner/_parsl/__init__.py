import logging
from pathlib import Path
from typing import Any
from typing import Dict
from typing import List

from parsl.dataflow.futures import AppFuture

from ...models import Workflow
from ..common import async_wrap
from ._setup import load_parsl_config


def get_app_future_result(app_future: AppFuture):
    """
    See issue #140 and https://stackoverflow.com/q/43241221/19085332

    By replacing
        .. = final_metadata.result()
    with
        .. = await async_wrap(get_app_future_result)(app_future=final_metadata)
    we avoid a (long) blocking statement.
    """
    return app_future.result()


def recursive_task_assembly(*args, **kwargs):
    raise NotImplementedError


async def process_workflow(
    *,
    workflow: Workflow,
    input_paths: List[Path],
    output_path: Path,
    input_metadata: Dict[str, Any],
    logger: logging.Logger,
    workflow_dir: Path,
    username: str = None,
) -> Dict[str, Any]:
    """
    Public interface to the runner backend

    Return
    ------
    final_metadata: Dict[str, Any]
        mapping representing the final state of the output dataset metadata
    """
    logger.info(f"{input_paths=}")
    logger.info(f"{output_path=}")

    # FIXME
    # in the following we most likely want a unique run id rather than the
    # generic name and id of the workflow
    with load_parsl_config(
        workflow_id=workflow.id,  # here
        workflow_name=workflow.name,  # here
        workflow_dir=workflow_dir,
        username=username,
        logger=logger,
    ) as dfk:
        final_metadata_future, dfk = recursive_task_assembly(
            data_flow_kernel=dfk,
            task=workflow,
            input_paths=input_paths,
            output_path=output_path,
            metadata=input_metadata,
            username=username,
        )
        logger.info("Definition of app futures complete, now start execution.")
        final_metadata = await async_wrap(get_app_future_result)(
            app_future=final_metadata_future
        )

    return final_metadata
