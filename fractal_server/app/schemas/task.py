from typing import Any
from typing import Dict
from typing import Optional

from sqlmodel import Field
from sqlmodel import SQLModel


__all__ = ("TaskCreate", "TaskUpdate", "TaskRead")


class _TaskBase(SQLModel):
    """
    Task base class

    A Task is the elemental unit of a workflow, and must be a self-standing
    executable.

    Attributes
    ----------
    name: str
        a human readable name for the task
    module: str
        module name of task
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
    module: str
    command: str
    source: str
    input_type: str
    output_type: str
    default_args: Dict[str, Any] = Field(default={})
    meta: Dict[str, Any] = Field(default={})

    class Config:
        arbitrary_types_allowed = True


class TaskUpdate(_TaskBase):
    name: Optional[str]  # type:ignore
    input_type: Optional[str]  # type:ignore
    output_type: Optional[str]  # type:ignore
    module: Optional[str]  # type:ignore
    command: Optional[str]  # type:ignore
    source: Optional[str]  # type:ignore
    default_args: Optional[Dict[str, Any]] = None  # type:ignore


class TaskCreate(_TaskBase):
    pass


class TaskRead(_TaskBase):
    id: int
