import logging
from concurrent.futures import Executor
from concurrent.futures import Future
from concurrent.futures import ThreadPoolExecutor
from functools import partial
from pathlib import Path
from typing import Any
from typing import Dict
from typing import List

from ...models import Workflow
from ...models import WorkflowTask
from .._common import call_single_parallel_task
from .._common import call_single_task
from ..common import TaskParameters


"""
Process Bakend

This backend runs fractal workflows as separate processes using a python
thread process pool, where each thread is responsible for running a single
task in a subprocess.

Incidentally, it represents the reference implementation for a backend.
"""


def call_parallel_task(
    *,
    executor: Executor,
    task: WorkflowTask,
    task_pars_depend_future: Future,  # py3.9 Future[TaskParameters],
    workflow_dir: Path,
) -> Future:  # py3.9 Future[TaskParameters]:
    """
    AKA collect results
    """
    task_pars_depend = task_pars_depend_future.result()
    component_list = task_pars_depend.metadata[task.parallelization_level]

    # Submit all tasks (one per component)
    partial_call_task = partial(
        call_single_parallel_task,
        task=task,
        task_pars=task_pars_depend,
        workflow_dir=workflow_dir,
    )
    map_iter = executor.map(partial_call_task, component_list)
    # Wait for execution of all parallel (this explicitly calls .result()
    # on each parallel task)
    for _ in map_iter:
        pass  # noqa: 701

    # Assemble a Future[TaskParameter]
    history = f"{task.task.name}: {component_list}"
    try:
        task_pars_depend.metadata["history"].append(history)
    except KeyError:
        task_pars_depend.metadata["history"] = [history]

    this_future: Future = Future()
    out_task_parameters = TaskParameters(
        input_paths=[task_pars_depend.output_path],
        output_path=task_pars_depend.output_path,
        metadata=task_pars_depend.metadata,
        logger=task_pars_depend.logger,
    )
    this_future.set_result(out_task_parameters)
    return this_future


def recursive_task_submission(
    *,
    executor: Executor,
    task_list: List[WorkflowTask],
    task_pars: TaskParameters,
    workflow_dir: Path,
) -> Future:
    """
    Recursively submit a list of task

    Each following task depends on the future.result() of the previous one,
    thus assuring the dependency chain.

    Induction process
    -----------------
    0: return a future which results in the task parameters necessary for the
       first task of the list

    n -> n+1: use output resulting from step `n` as task parameters to submit
       task `n+1`

    Return
    ------
    this_future (Future[TaskParameters]):
        a future that results to the task parameters which constitute the
        input of the following task in the list.
    """
    try:
        *dependencies, this_task = task_list
    except ValueError:
        # step 0: return future containing original task_pars
        pseudo_future: Future = Future()
        pseudo_future.set_result(task_pars)
        return pseudo_future

    # step n => step n+1
    task_pars.logger.debug(f"submitting task {this_task.order=}")
    # parallelization_level = this_task.task.parallelization_level

    task_pars_depend_future = recursive_task_submission(
        executor=executor,
        task_list=dependencies,
        task_pars=task_pars,
        workflow_dir=workflow_dir,
    )

    if this_task.is_parallel:
        this_future = call_parallel_task(
            executor=executor,
            task=this_task,
            task_pars_depend_future=task_pars_depend_future,
            workflow_dir=workflow_dir,
        )
    else:
        this_future = executor.submit(
            call_single_task,
            task=this_task,
            task_pars=task_pars_depend_future.result(),
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
    TODO:
    in case of failure we must return the most recent clean metadata

    Returns:
    output_dataset_metadata (Dict):
        the output metadata
    """

    with ThreadPoolExecutor() as executor:
        output_dataset_metadata = recursive_task_submission(
            executor=executor,
            task_list=workflow.task_list,
            task_pars=TaskParameters(
                input_paths=input_paths,
                output_path=output_path,
                metadata=input_metadata,
                logger=logger,
            ),
            workflow_dir=workflow_dir,
        )
    return output_dataset_metadata.result()
