import logging
from concurrent.futures import Future
from pathlib import Path
from typing import Any
from typing import Dict
from typing import List

from parsl.app.app import join_app
from parsl.app.app import python_app
from parsl.app.python import PythonApp
from parsl.dataflow.dflow import DataFlowKernel
from parsl.dataflow.futures import AppFuture

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
) -> AppFuture:  # AppFuture[TaskParameters]

    # TODO: add executor
    _this_parallel_component_app = PythonApp(
        _this_parallel_component, data_flow_kernel=data_flow_kernel
    )

    # TODO: add executor?
    _collect_results_and_assemble_history_app = PythonApp(
        _collect_results_and_assemble_history,
        data_flow_kernel=data_flow_kernel,
    )

    # TODO: add executor?
    @join_app(data_flow_kernel=data_flow_kernel)
    def _parallel_task_app_future_app(task_pars) -> AppFuture:
        component_list = task_pars.metadata[task.parallelization_level]
        # app che colleziona i risultati
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
) -> AppFuture:  # AppFuture[TaskParameters]
    if not workflow_dir:
        raise RuntimeError

    # TODO: add executor?
    # assemble full args
    @python_app(data_flow_kernel=data_flow_kernel)
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

    # TODO: use executor

    try:
        *dependencies, this_task = task_list
    except ValueError:
        # step 0: return future containing original task_pars
        pseudo_future: Future = Future()
        pseudo_future.set_result(task_pars)
        return pseudo_future
    # step n => step n+1
    logger.debug(f"submitting task {this_task.order=}")
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
            task_pars_depend_future=task_pars_depend_future,
            workflow_dir=workflow_dir,
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
        logger.info("Definition of app futures complete, now start execution.")
        final_metadata = await async_wrap(get_app_future_result)(
            app_future=final_metadata_future
        )

    return final_metadata
