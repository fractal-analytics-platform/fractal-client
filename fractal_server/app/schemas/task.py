from typing import Any
from typing import Dict
from typing import List

from sqlmodel import Field
from sqlmodel import SQLModel

from ..schemas import TaskRead


class WorkflowBase(SQLModel):
    name: str
    project_id: int


class WorkflowRead(WorkflowBase):
    id: int
    task_list: List[TaskRead]


class TaskBase(SQLModel):
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

    class Config:
        arbitrary_types_allowed = True
