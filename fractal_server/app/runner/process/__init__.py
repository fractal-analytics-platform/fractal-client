import json
import logging
import subprocess  # nosec
from concurrent.futures import Executor
from concurrent.futures import Future
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from shlex import split as shlex_split
from typing import Any
from typing import Dict
from typing import List

from ...models import Workflow
from ...models import WorkflowTask
from ..common import TaskParameterEncoder
from ..common import TaskParameters


"""
Process Bakend

This backend runs fractal workflows as separate processes using a python
process pool.

Incidentally, it represents the reference implementation for a backend.
"""


class TaskExecutionError(RuntimeError):
    def __init__(self, completed_process: subprocess.CompletedProcess):
        self.completed_process = completed_process
        super().__init__(completed_process.stderr.decode("utf-8"))


def _call_command_wrapper(cmd: str) -> subprocess.CompletedProcess:
    """
    Call command and return stdout, stderr, retcode
    """

    result = subprocess.run(shlex_split(cmd), capture_output=True)  # nosec
    if result.returncode != 0:
        raise TaskExecutionError(result)
    return result


def call_single_task(
    task: WorkflowTask,
    task_pars: TaskParameters,
    workflow_dir: Path = None,
) -> TaskParameters:
    from devtools import debug

    if not workflow_dir:
        raise RuntimeError
    # assemble full args
    args_dict = task_pars.dict(exclude={"logger"})
    debug(task.arguments)
    args_dict.update(task.arguments)
    debug(args_dict)

    # write args file
    args_file_path = workflow_dir / f"{task.order}.args.json"
    with open(args_file_path, "w") as f:
        json.dump(args_dict, f, cls=TaskParameterEncoder)

    # assemble full command
    debug(task.task)
    cmd = f"{task.task.command} -j {args_file_path}"
    debug(cmd)
    return _call_command_wrapper(cmd)


def recursive_task_submission(
    executor: Executor,
    task_list: List[WorkflowTask],
    task_pars: TaskParameters,
) -> Future:
    """
    0: return a future with everything that is needed to submit first task of
       workflow

    n -> n+1: use output from n to submit task n+1
    """
    from devtools import debug

    debug(f"TASK SUBMISSION: {[t.order for t in task_list]}")
    try:
        *dependent_tasks, this_task = task_list
    except ValueError:
        # 0
        debug("0")
        pseudo_future = Future()
        pseudo_future.set_result(task_pars)
        return pseudo_future

    # n -> n+1
    return executor.submit(
        call_single_task,
        task=this_task,
        task_pars=recursive_task_submission(
            executor=executor, task_list=dependent_tasks
        ).result(),
    )


async def process_workflow(
    *,
    workflow: Workflow,
    input_paths: List[Path],
    output_path: Path,
    input_metadata: Dict[str, Any],
    logger: logging.Logger,
    username: str = None,
) -> Dict[str, Any]:

    from devtools import debug

    debug(workflow)

    with ProcessPoolExecutor() as executor:
        recursive_task_submission(
            executor=executor,
            task_list=workflow.task_list,
            task_pars=TaskParameters(
                input_paths=input_paths,
                output_path=output_path,
                metadata=input_metadata,
                logger=logger,
            ),
        )

    raise NotImplementedError
