from typing import Any
from typing import Dict
from typing import List
from typing import Optional

from sqlalchemy import Column
from sqlalchemy.types import JSON
from sqlmodel import Field
from sqlmodel import Relationship
from sqlmodel import SQLModel

from ..schemas.task import WorkflowBase
from .task import Task


class LinkTaskWorkflow(SQLModel, table=True):
    """
    Crossing table between Task and Workflow

    In addition to the foreign keys, it allows for parameter overriding and
    keeps the order within the list of tasks of the workflow.

    Attributes
    ----------
    """

    class Config:
        arbitrary_types_allowed = True
        fields = {"parent": {"exclude": True}}

    id: Optional[int] = Field(default=None, primary_key=True)

    workflow_id: Optional[int] = Field(foreign_key="workflow.id")
    task_id: Optional[int] = Field(foreign_key="task.id")

    order: Optional[int]
    args: Dict[str, Any] = Field(sa_column=Column(JSON), default={})


class Workflow(WorkflowBase, table=True):
    """
    Workflow

    Attributes
    ----------
    task_list: List(LinkTaskWorkflow)
        List of associations to tasks.
    """

    id: Optional[int] = Field(default=None, primary_key=True)
    project_id: int = Field(foreign_key="project.id")

    task_list: List["LinkTaskWorkflow"] = Relationship(
        sa_relationship_kwargs=dict(
            lazy="selectin",
            order_by="LinkTaskWorkflow.order",
            collection_class=ordering_list("order"),
        ),
    )

    def insert_task(
        self,
        task: Task,
        *,
        args: Dict[str, Any] = None,
        order: Optional[int] = None
    ) -> None:
        if order is None:
            order = len(self.task_list)

        self.task_list.insert(order, LinkTaskWorkflow(task_id=task.id))
