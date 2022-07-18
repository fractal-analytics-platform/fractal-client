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
from ..models.project import Resource
from ..models.task import PreprocessedTask
from ..models.task import Subtask
from ..models.task import Task


@parsl.python_app
def _collect_results(inputs: List[PythonApp]):
    return [x for x in inputs]


parsl.load(Config())


def _process_workflow(
    task: Union[Task, Subtask],
    input_paths: List[Path],
    output_path: Path,
    intermediate_path: Optional[Path] = None,
) -> PythonApp:
    """
    Unwrap workflow recursively

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

    apps = []

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


async def submit_workflow(
    *,
    input_dataset: Dataset,
    workflow: Task,
    output_dataset: Optional[Union[Dataset, str]] = None,
    output_resource: Optional[Union[Resource, str]] = None,
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

    if workflow.input_type != input_dataset.type:
        raise TypeError(
            f"Incompatible types `{workflow.input_type}` of workflow "
            f"`{workflow.name}` and `{input_dataset.type}` of dataset "
            f"`{input_dataset.name}`"
        )

    input_path = [r.glob_path for r in input_dataset.resource_list]

    output_path = None
    metadata = input_dataset.metadata

    if "workflow" in workflow.resource_type:
        for subtask in workflow.subtask_list:
            kwargs = subtask._arguments

            @parsl.python_app()
            def workflow_app(
                input_path: List[Path] = input_path,
                output_path: Path = output_path,
                metadata: Dict[str, Any] = metadata,
                kwargs: Dict[str, Any] = kwargs,
            ):
                task_module = importlib.import_module(subtask.import_path)
                __callable = getattr(task_module, subtask.callable)
                __callable(
                    input_path=input_path,
                    output_path=output_dataset,
                    metadata=metadata,
                    **kwargs,
                )

    app_list = []
    # TODO
    # for parallelization_item in parallelization_level:
    #     app_list.append(workflow_app(parallelization_item))

    map(lambda app: app.result(), app_list)
