import json
import logging
import subprocess  # nosec
from pathlib import Path
from shlex import split as shlex_split

from ..models import WorkflowTask
from .common import TaskParameterEncoder
from .common import TaskParameters
from .common import write_args_file


class TaskExecutionError(RuntimeError):
    pass


def _call_command_wrapper(cmd: str, stdout: Path, stderr: Path) -> None:
    """
    Call command and return stdout, stderr, retcode
    """
    fp_stdout = stdout.open("w")
    fp_stderr = stderr.open("w")
    try:
        result = subprocess.run(  # nosec
            shlex_split(cmd),
            stderr=fp_stderr,
            stdout=fp_stdout,
        )
    except Exception as e:
        raise e
    finally:
        fp_stdout.close()
        fp_stderr.close()
    if result.returncode != 0:
        with stderr.open("r") as fp_stderr:
            err = fp_stderr.read()
        raise TaskExecutionError(err)


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

    logger = logging.getLogger(task_pars.logger_name)

    stdout_file = workflow_dir / f"{task.order}.out"
    stderr_file = workflow_dir / f"{task.order}.err"
    metadata_diff_file = workflow_dir / f"{task.order}.metadiff.json"

    # assemble full args
    args_dict = task.assemble_args(extra=task_pars.dict())

    # write args file
    args_file_path = workflow_dir / f"{task.order}.args.json"
    write_args_file(args=args_dict, path=args_file_path)

    # assemble full command
    cmd = (
        f"{task.task.command} -j {args_file_path} "
        f"--metadata-out {metadata_diff_file}"
    )

    logger.debug(f"executing task {task.order=}")
    _call_command_wrapper(cmd, stdout=stdout_file, stderr=stderr_file)

    # NOTE:
    # This assumes that the new metadata is printed to stdout
    # and nothing else outputs to stdout
    with metadata_diff_file.open("r") as f_metadiff:
        diff_metadata = json.load(f_metadiff)
    updated_metadata = task_pars.metadata.copy()
    updated_metadata.update(diff_metadata)

    # Assemble a Future[TaskParameter]
    history = f"{task.task.name}"
    try:
        updated_metadata["history"].append(history)
    except KeyError:
        updated_metadata["history"] = [history]

    out_task_parameters = TaskParameters(
        input_paths=[task_pars.output_path],
        output_path=task_pars.output_path,
        metadata=updated_metadata,
        logger_name=task_pars.logger_name,
    )
    return out_task_parameters


def call_single_parallel_task(
    component: str,
    *,
    task: WorkflowTask,
    task_pars: TaskParameters,
    workflow_dir: Path = None,
) -> None:
    if not workflow_dir:
        raise RuntimeError
    logger = logging.getLogger(task_pars.logger_name)

    prefix = f"{task.order}_par_{component}"
    stdout_file = workflow_dir / f"{prefix}.out"
    stderr_file = workflow_dir / f"{prefix}.err"
    metadata_diff_file = workflow_dir / f"{prefix}.metadiff.json"

    logger.debug(f"calling task {task.order=} on {component=}")
    # FIXME refactor with `write_args_file` and `task.assemble_args`
    # assemble full args
    args_dict = task_pars.dict()
    args_dict.update(task.arguments)
    args_dict["component"] = component

    # write args file
    args_file_path = workflow_dir / f"{prefix}.args.json"
    with open(args_file_path, "w") as f:
        json.dump(args_dict, f, cls=TaskParameterEncoder)
    # FIXME: UP TO HERE

    # assemble full command
    cmd = (
        f"{task.task.command} -j {args_file_path} "
        f"--metadata-out {metadata_diff_file}"
    )

    logger.debug(f"executing task {task.order=}")
    _call_command_wrapper(cmd, stdout=stdout_file, stderr=stderr_file)
