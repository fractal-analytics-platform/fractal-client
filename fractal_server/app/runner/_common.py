import json
import subprocess  # nosec
from pathlib import Path
from shlex import split as shlex_split

from ..models import WorkflowTask
from .common import TaskParameters
from .common import write_args_file


class TaskExecutionError(RuntimeError):
    """
    Indicate that the subprocess execution exited status != 0

    The traceback is extracted from the subprocess stderr and used to
    initialise the error. If the command is a Python executable, this gives
    access to the full traceback.

    Attributes
    ----------
    completed_process (subprocess.CompletedProcess):
        the full object as returned by subprocess.run()
    """

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
    *,
    task: WorkflowTask,
    task_pars: TaskParameters,
    workflow_dir: Path = None,
) -> TaskParameters:
    """
    Call a single task

    This assemble the runner (input_paths, output_path, ...) and task
    arguments (arguments that are specific to the task, such as message or
    index in the dummy task), writes them to file, call the task executable
    command passing the arguments file as an input and assembles the output

    Attributes
    ----------
    task (WorkflowTask):
        the workflow task to be called. This includes task specific arguments
        via the task.task.arguments attribute.
    task_pars (TaskParameters):
        the parameters required to run the task which are not specific to the
        task, e.g., I/O paths.
    workflow_dir (Path):
        the directory in which the execution takes place, and where all
        artifacts are written.

    Return
    ------
    out_task_parameters (TaskParameters):
        a TaskParameters in which the previous output becomes the input and
        where metadata is the metadata dictionary returned by the task being
        called.
    """
    if not workflow_dir:
        raise RuntimeError

    # assemble full args
    args_dict = task.assemble_args(extra=task_pars.dict(exclude={"logger"}))

    # write args file
    args_file_path = workflow_dir / f"{task.order}.args.json"
    write_args_file(args=args_dict, path=args_file_path)

    # assemble full command
    cmd = f"{task.task.command} -j {args_file_path}"

    task_pars.logger.debug(f"executing task {task.order=}")
    completed_process = _call_command_wrapper(cmd)

    # NOTE:
    # This assumes that the new metadata is printed to stdout
    # and nothing else outputs to stdout
    diff_metadata = json.loads(completed_process.stdout)
    updated_metadata = task_pars.metadata.copy()
    updated_metadata.update(diff_metadata)

    out_task_parameters = TaskParameters(
        input_paths=[task_pars.output_path],
        output_path=task_pars.output_path,
        metadata=updated_metadata,
        logger=task_pars.logger,
    )
    return out_task_parameters
