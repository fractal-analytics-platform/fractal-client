import logging
from concurrent.futures import Future
from pathlib import Path
from typing import Any
from typing import Dict
from typing import List
from typing import Literal
from typing import Union

from parsl.app.app import join_app
from parsl.app.app import python_app
from parsl.app.python import PythonApp
from parsl.dataflow.dflow import DataFlowKernel
from parsl.dataflow.futures import AppFuture

from ....config_runner import settings
from ...models import Workflow
from ...models import WorkflowTask
from .._common import call_single_task
from ..common import async_wrap
from ..common import TaskParameters
from ._functions import _collect_results_and_assemble_history
from ._functions import _this_parallel_component
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
    executors: Union[List[str], Literal["all"]],
) -> AppFuture:  # AppFuture[TaskParameters]

    # Define PythonApp for execution of sigle instance of parallel task
    _this_parallel_component_app = PythonApp(
        _this_parallel_component,
        data_flow_kernel=data_flow_kernel,
        executors=executors,
    )

    # Define PythonApp for execution of sigle instance of parallel task
    _collect_results_and_assemble_history_app = PythonApp(
        _collect_results_and_assemble_history,
        data_flow_kernel=data_flow_kernel,
        executors=executors,
    )

    @join_app(data_flow_kernel=data_flow_kernel)
    def _parallel_task_app_future_app(task_pars) -> AppFuture:
        component_list = task_pars.metadata[task.parallelization_level]
        dependency_futures = [
            _this_parallel_component_app(
                task_pars=task_pars,
                component=c,
                task=task,
                workflow_dir=workflow_dir,
            )
            for c in component_list
        ]

        return _collect_results_and_assemble_history_app(
            task_pars=task_pars,
            component_list=component_list,
            inputs=dependency_futures,
            task=task,
        )

    return _parallel_task_app_future_app(task_pars_depend_future)


def _serial_task_assembly(
    data_flow_kernel: DataFlowKernel,
    task: WorkflowTask,
    task_pars_depend_future: AppFuture,
    workflow_dir: Path,
    executors: Union[List[str], Literal["all"]],
) -> AppFuture:  # AppFuture[TaskParameters]
    if not workflow_dir:
        raise RuntimeError

    # assemble full args
    @python_app(data_flow_kernel=data_flow_kernel, executors=executors)
    def _this_app(task_pars):
        return call_single_task(
            task=task,
            task_pars=task_pars,
            workflow_dir=workflow_dir,
        )

    app_future = _this_app(task_pars=task_pars_depend_future)
    return app_future


def recursive_task_assembly(
    *,
    data_flow_kernel: DataFlowKernel,
    task_list: List[WorkflowTask],
    task_pars: TaskParameters,
    workflow_dir: Path,
) -> AppFuture:

    logger = logging.getLogger(task_pars.logger_name)

    try:
        *dependencies, this_task = task_list
    except ValueError:
        # step 0: return future containing original task_pars
        pseudo_future: Future = Future()
        pseudo_future.set_result(task_pars)
        return pseudo_future

    # Extract and validate task executor
    executors: Union[List[str], Literal["all"]]
    if this_task.executor:
        # Verify match between new_task_executor and available executors
        if this_task.executor not in data_flow_kernel.executors.keys():
            msg = (
                f"{this_task.executor=} is not in "
                f"{data_flow_kernel.executors.keys()=}"
            )
            logger.error(msg)
            raise ValueError(msg)
        executors = [this_task.executor] if this_task.executor else "all"
    else:
        executors = [settings.RUNNER_DEFAULT_EXECUTOR]

    # step n => step n+1
    logger.debug(f"submitting task {this_task.order=} to {executors=}")
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
            executors=executors,
        )
    else:
        this_future = _serial_task_assembly(
            data_flow_kernel=data_flow_kernel,
            task=this_task,
            task_pars_depend_future=task_pars_depend_future,
            workflow_dir=workflow_dir,
            executors=executors,
        )

    return this_future


async def process_workflow(
    *,
    workflow: Workflow,
    input_paths: List[Path],
    output_path: Path,
    input_metadata: Dict[str, Any],
    logger_name: str,
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
    logger = logging.getLogger(logger_name)
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
        logger_name=logger_name,
    ) as dfk:
        logger.info("Start definition of app futures (from last to first)")
        final_metadata_future = recursive_task_assembly(
            data_flow_kernel=dfk,
            task_list=workflow.task_list,
            task_pars=TaskParameters(
                input_paths=input_paths,
                output_path=output_path,
                metadata=input_metadata,
                logger_name=logger_name,
            ),
            workflow_dir=workflow_dir,
        )
        logger.info("Definition of app futures complete")
        logger.info("Start execution")
        final_metadata = await async_wrap(get_app_future_result)(
            app_future=final_metadata_future
        )
        logger.info("Execution complete")

    return final_metadata
