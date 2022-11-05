from pathlib import Path
from typing import Any
from typing import Dict
from typing import List
from typing import Literal
from typing import Optional

from pydantic import BaseModel
from sqlmodel import Field  # type: ignore
from sqlmodel import SQLModel


__all__ = (
    "TaskCreate",
    "TaskUpdate",
    "TaskRead",
    "TaskCollectPip",
    "TaskCollectResult",
    "TaskCollectStatus",
)


class _TaskBase(SQLModel):
    """
    Task base class

    A Task is the elemental unit of a workflow, and must be a self-standing
    executable.

    Attributes
    ----------
    name: str
        a human readable name for the task
    command: str
        the command(s) that executes the task
    source: str
        path or url to task source. This is the information is used to match
        tasks across fractal installations when a workflow is imported.
    input_type, output_type: str
        the type of data the task expects as input, output, respectively.
    default_args: Dict[str, Any]
        dictionary (saved as JSON) of the default parameters of the task
    """

    name: str
    command: str
    source: str
    input_type: str
    output_type: str
    default_args: Dict[str, Any] = Field(default={})
    meta: Optional[Dict[str, Any]] = Field(default={})

    class Config:
        arbitrary_types_allowed = True


class TaskUpdate(_TaskBase):
    name: Optional[str]  # type:ignore
    input_type: Optional[str]  # type:ignore
    output_type: Optional[str]  # type:ignore
    command: Optional[str]  # type:ignore
    source: Optional[str]  # type:ignore
    default_args: Optional[Dict[str, Any]] = None  # type:ignore


class TaskRead(_TaskBase):
    id: int


class TaskCreate(_TaskBase):
    pass


class _TaskCollectBase(BaseModel):
    pass


class TaskCollectPip(_TaskCollectBase):
    package: str
    version: Optional[str]
    python_version: str = "3.8"
    package_extras: Optional[str]


class TaskCollectResult(_TaskCollectBase):
    """
    Temporary result of the collection process

    This schema represents the partial result returned when collection is
    initiated but it is continuing in the background.

    Attributes:
        package (str):
            the package name
        venv_path (Path):
            the installation path relative to `FRACTAL_ROOT`
        info (str):
            arbitrary information
    """

    package: str
    venv_path: Path
    info: Optional[str]


class TaskCollectStatus(_TaskCollectBase):
    status: Literal["OK", "pending", "fail"]
    task_list: Optional[List[TaskRead]] = Field(default=[])
    log: Optional[str]
