import asyncio
import json
from functools import partial
from functools import wraps
from json import JSONEncoder
from pathlib import Path
from typing import Any
from typing import Callable
from typing import Dict
from typing import List
from typing import Optional

from pydantic import BaseModel

from ...utils import close_logger as close_job_logger  # noqa F401
from ..models import Dataset
from ..models import Project
from ..models.task import Task


class TaskParameterEncoder(JSONEncoder):
    def default(self, value):
        if isinstance(value, Path):
            return value.as_posix()
        return JSONEncoder.default(self, value)


class TaskParameters(BaseModel):
    input_paths: List[Path]
    output_path: Path
    metadata: Dict[str, Any]
    logger_name: Optional[str] = None
    username: Optional[str] = None

    class Config:
        arbitrary_types_allowed = True
        extra = "forbid"


async def auto_output_dataset(
    *,
    project: Project,
    input_dataset: Dataset,
    workflow: Task,
    overwrite_input: bool = False,
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

        if len(output_dataset.paths) != 1:
            raise ValueError(
                "Cannot determine output path: Multiple paths in dataset."
            )
        else:
            output_path = output_dataset.paths[0]
    return output_path


def async_wrap(func: Callable) -> Callable:
    """
    See issue #140 and https://stackoverflow.com/q/43241221/19085332

    By replacing
        .. = final_metadata.result()
    with
        .. = await async_wrap(get_app_future_result)(app_future=final_metadata)
    we avoid a (long) blocking statement.
    """

    @wraps(func)
    async def run(*args, loop=None, executor=None, **kwargs):
        if loop is None:
            loop = asyncio.get_event_loop()
        pfunc = partial(func, *args, **kwargs)
        return await loop.run_in_executor(executor, pfunc)

    return run


def write_args_file(args: Dict[str, Any], path: Path):
    with path.open("w") as f:
        json.dump(args, f, cls=TaskParameterEncoder, indent=4)
