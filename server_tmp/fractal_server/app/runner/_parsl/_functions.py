from pathlib import Path
from typing import List

from parsl.dataflow.futures import AppFuture

from ...models import WorkflowTask
from .._common import call_single_parallel_task
from ..common import TaskParameters


def _this_parallel_component(
    task_pars: TaskParameters,
    component: str,
    task: WorkflowTask,
    workflow_dir: Path,
) -> None:
    return call_single_parallel_task(
        component,
        task=task,
        task_pars=task_pars,
        workflow_dir=workflow_dir,
    )


def _collect_results_and_assemble_history(
    task_pars: TaskParameters,
    component_list: List[str],
    task: WorkflowTask,
    inputs: List[AppFuture],
) -> TaskParameters:
    history = f"{task.task.name}: {component_list}"
    try:
        task_pars.metadata["history"].append(history)
    except KeyError:
        task_pars.metadata["history"] = [history]

    out_task_parameters = TaskParameters(
        input_paths=[task_pars.output_path],
        output_path=task_pars.output_path,
        metadata=task_pars.metadata,
        logger_name=task_pars.logger_name,
    )
    return out_task_parameters
