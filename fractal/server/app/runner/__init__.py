import importlib
from pathlib import Path
from typing import Any
from typing import Dict
from typing import List
from typing import Optional
from typing import Union

import parsl
from parsl.app.python import PythonApp
from parsl.config import Config
from parsl.dataflow.futures import AppFuture

from ..models.project import Dataset
from ..models.project import Project
from ..models.task import PreprocessedTask
from ..models.task import Subtask
from ..models.task import Task


@parsl.python_app
def _collect_results(inputs: List[PythonApp]):
    return [x for x in inputs]


parsl.load(Config())


def _atomic_task_factory(
    *,
    task: Union[Task, Subtask, PreprocessedTask],
    input_paths: List[Path],
    output_path: Path,
    metadata: Optional[Dict[str, Any]] = None,
    depends_on: Optional[List[AppFuture]] = None,
) -> PythonApp:
    """
    Single task processing

    Create a single PARSL app that encapsulates the task at hand and
    its parallelizazion.
    """
    if depends_on is None:
        depends_on = []

    task_args = task._arguments

    # TODO
    # Pass executor
    # executor = task_args.get("executor", "cpu")
    # @parsl.python_app(executors=[executor])

    @parsl.python_app
    def _task_app(
        *,
        input_paths: List[Path] = input_paths,
        output_path: Path = output_path,
        metadata: Optional[Dict[str, Any]] = metadata,
        task_args: Optional[Dict[str, Any]] = task_args,
        component: Optional[Dict[str, Any]] = None,
        inputs=None,
    ):
        if component is None:
            component = {}

        task_module = importlib.import_module(task.import_path)
        _callable = getattr(task_module, task.callable)
        return _callable(
            input_paths=input_paths,
            output_path=output_path,
            metadata=metadata,
            **component,
            **task_args,
        )

    parall_level = task_args.get("parallelization_level", None)
    if metadata and parall_level:
        parall_item_gen = (par_item for par_item in metadata[parall_level])
        return _collect_results(
            inputs=[
                _task_app(component={parall_level: item})
                for item in parall_item_gen
            ]
        )
    else:
        return _task_app(inputs=depends_on)


def _process_workflow(
    task: Union[Task, Subtask],
    input_paths: List[Path],
    output_path: Path,
) -> PythonApp:
    """
    Creates the PARSL app that will execute the full workflow, taking care of
    dependencies

    Arguments
    ---------
    output_path (Path):
        directory or file where the final output, i.e., the output of the last
        task, will be written
    TBD

    Return
    ------
    TBD
    """
    preprocessed = task.preprocess()

    this_input = input_paths
    this_output = output_path

    apps: List[PythonApp] = []

    for i, task in enumerate(preprocessed):
        this_task_app = _atomic_task_factory(
            task=task,
            input_paths=this_input,
            output_path=this_output,
            depends_on=[apps[i - 1] if i > 0 else None],
        )
        apps.append(this_task_app)
        this_input = [this_output]

    # Got to make sure that it is executed serially, task by task
    return apps[-1]


async def auto_output_dataset(
    *,
    project: Project,
    input_dataset: Dataset,
    workflow: Task,
    overwrite_input: bool = False
):
    """
    Determine the output dataset if it was not provided explicitly

    Only datasets containing exactly one path can be used as output.

    Returns
    -------
    output_dataset (Dataset):
        the output dataset
    """
    if overwrite_input and not input_dataset.read_only:
        input_paths = input_dataset.paths
        if len(input_paths) != 1:
            raise ValueError
        output_dataset = input_dataset
    else:
        raise NotImplementedError

    return output_dataset


def validate_workflow_compatibility(
    *,
    input_dataset: Dataset,
    workflow: Task,
    output_dataset: Optional[Dataset] = None,
):
    """
    Check compatibility of workflow and input / ouptut dataset
    """
    if (
        workflow.input_type != "Any"
        and workflow.input_type != input_dataset.type
    ):
        raise TypeError(
            f"Incompatible types `{workflow.input_type}` of workflow "
            f"`{workflow.name}` and `{input_dataset.type}` of dataset "
            f"`{input_dataset.name}`"
        )

    if not output_dataset:
        if input_dataset.read_only:
            raise ValueError("Input dataset is read-only")
        else:
            input_paths = input_dataset.paths
            if len(input_paths) != 1:
                # Only single input can be safely transformed in an output
                raise ValueError(
                    "Cannot determine output path: multiple input "
                    "paths to overwrite"
                )
            else:
                output_path = input_paths[0]
    else:
        output_path = output_dataset.paths
        if len(output_path) != 1:
            raise ValueError(
                "Cannot determine output path: Multiple paths in dataset."
            )
    return output_path


async def submit_workflow(
    *,
    workflow: Task,
    input_dataset: Dataset,
    output_dataset: Dataset,
):
    """
    Arguments
    ---------
    output_dataset (Dataset | str) :
        the destination dataset of the workflow. If not provided, overwriting
        of the input dataset is implied and an error is raised if the dataset
        is in read only mode. If a string is passed and the dataset does not
        exist, a new dataset with that name is created and within it a new
        resource with the same name.
    """
    if output_dataset is None:
        if input_dataset.read_only:
            raise ValueError(
                "Input dataset is read-only: cannot overwrite. "
                "Please provide `output_dataset` explicitly."
            )
        else:
            output_dataset = input_dataset
    else:
        pass
        # TODO: check that dataset exists and if not create dataset and
        # resource.

    input_paths = input_dataset.paths
    output_path = output_dataset.paths[0]

    wf = _process_workflow(
        task=workflow,
        input_paths=input_paths,
        output_path=output_path,
    )
    from devtools import debug
    debug(wf)
    from time import sleep
    sleep(2)
    debug(wf)
    raise NotImplementedError
