from typing import Any
from typing import Dict
from typing import List
from typing import Optional

from sqlalchemy import Column
from sqlalchemy.ext.orderinglist import ordering_list
from sqlalchemy.types import JSON
from sqlmodel import Field
from sqlmodel import Relationship
from sqlmodel import SQLModel

from ...common.models.task import WorkflowBase
from .task import Task


class LinkTaskWorkflow(SQLModel, table=True):
    """
    Crossing table between Task and Workflow

    In addition to the foreign keys, it allows for parameter overriding.

    Attributes
    ----------
    """

    class Config:
        arbitrary_types_allowed = True
        fields = {"parent": {"exclude": True}}

    id: Optional[int] = Field(default=None, primary_key=True)

    task_id: Optional[int] = Field(default=None, foreign_key="task.id")
    task: "Task" = Relationship()

    order: Optional[int]
    args: Dict[str, Any] = Field(sa_column=Column(JSON), default={})


class Workflow(WorkflowBase, table=True):
    id: Optional[int] = Field(default=None, primary_key=True)
    project_id: int = Field(foreign_key="project.id")

