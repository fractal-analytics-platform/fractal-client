import logging
from concurrent.futures import Future
from pathlib import Path
from typing import Any
from typing import Dict
from typing import List

from parsl.app.app import bash_app
from parsl.dataflow.dflow import DataFlowKernel
from parsl.dataflow.futures import AppFuture

from ...models import Workflow
from ...models import WorkflowTask
from ..common import async_wrap
from ..common import TaskParameters
from ..common import write_args_file
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


def _parallel_task_assembly(
    data_flow_kernel: DataFlowKernel,
    task: WorkflowTask,
    task_pars_depend_future: AppFuture,
    workflow_dir: Path,
    parallelization_level: str,
) -> AppFuture:
    raise NotImplementedError


def _serial_task_assembly(
    data_flow_kernel: DataFlowKernel,
    task: WorkflowTask,
    task_pars_depend_future: AppFuture,
    workflow_dir: Path,
) -> AppFuture:
    if not workflow_dir:
        raise RuntimeError

    stdout_file = workflow_dir / f"{task.order}.out.json"
    stderr_file = workflow_dir / f"{task.order}.err"

    # assemble full args
    @bash_app(data_flow_kernel=data_flow_kernel)
    def _this_bash_app(inputs, stderr=stderr_file.as_posix()):

        task_pars = inputs[0]

        args_dict = task.assemble_args(
            extra=task_pars.dict(exclude={"logger"})
        )

        # write args file
        args_file_path = workflow_dir / f"{task.order}.args.json"
        write_args_file(args=args_dict, path=args_file_path)

        return f"{task.task.command} -j {args_file_path} -o {stdout_file}"

    app_future = _this_bash_app(inputs=[task_pars_depend_future])
    return app_future


def recursive_task_assembly(
    *,
    data_flow_kernel: DataFlowKernel,
    task_list: List[WorkflowTask],
    task_pars: TaskParameters,
    workflow_dir: Path,
) -> AppFuture:

    try:
        *dependencies, this_task = task_list
    except ValueError:
        # step 0: return future containing original task_pars
        pseudo_future: Future = Future()
        pseudo_future.set_result(task_pars)
        return pseudo_future
    # step n => step n+1
    task_pars.logger.debug(f"submitting task {this_task.order=}")
    parallelization_level = this_task.task.parallelization_level

    task_pars_depend_future = recursive_task_assembly(
        data_flow_kernel=data_flow_kernel,
        task_list=dependencies,
        task_pars=task_pars,
        workflow_dir=workflow_dir,
    )

    if parallelization_level:
        this_future = _parallel_task_assembly(
            data_flow_kernel=data_flow_kernel,
            task=this_task,
            task_pars_depend_future=task_pars_depend_future,
            workflow_dir=workflow_dir,
            parallelization_level=parallelization_level,
        )
    else:
        this_future = _serial_task_assembly(
            data_flow_kernel=data_flow_kernel,
            task=this_task,
            task_pars_depend_future=task_pars_depend_future.result(),
            workflow_dir=workflow_dir,
        )
    return this_future


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
        final_metadata_future = recursive_task_assembly(
            data_flow_kernel=dfk,
            task_list=workflow.task_list,
            task_pars=TaskParameters(
                input_paths=input_paths,
                output_path=output_path,
                metadata=input_metadata,
                logger=logger,
            ),
            workflow_dir=workflow_dir,
        )
        logger.info("Definition of app futures complete, now start execution.")
        final_metadata = await async_wrap(get_app_future_result)(
            app_future=final_metadata_future
        )

    return final_metadata
