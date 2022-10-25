from typing import Any
from typing import Dict
from typing import List
from typing import Optional

from sqlalchemy import Column
from sqlalchemy.ext.orderinglist import ordering_list
from sqlalchemy.types import JSON
from sqlmodel import Field
from sqlmodel import Relationship

from ..db import AsyncSession
from ..schemas.workflow import _WorkflowBase
from ..schemas.workflow import _WorkflowTaskBase
from .models_utils import popget
from .task import Task


class WorkflowTask(_WorkflowTaskBase, table=True):
    """
    A Task as part of a Workflow

    This is a crossing table between Task and Workflow. In addition to the
    foreign keys, it allows for parameter overriding and keeps the order
    within the list of tasks of the workflow.

    Attributes
    ----------
    TODO
    """

    class Config:
        arbitrary_types_allowed = True
        fields = {"parent": {"exclude": True}}

    id: Optional[int] = Field(default=None, primary_key=True)

    workflow_id: Optional[int] = Field(foreign_key="workflow.id")
    task_id: Optional[int] = Field(foreign_key="task.id")

    order: Optional[int]
    args: Dict[str, Any] = Field(sa_column=Column(JSON), default={})

    task: Task = Relationship(sa_relationship_kwargs=dict(lazy="selectin"))

    @property
    def arguments(self):
        """
        Override default arguments and strip specific arguments (executor and
        parallelization_level)
        """
        out = self.task.default_args.copy()
        out.update(self.args)
        popget(out, "parallelization_level")
        return out


class Workflow(_WorkflowBase, table=True):
    """
    Workflow

    Attributes
    ----------
    task_list: List(LinkTaskWorkflow)
        List of associations to tasks.
    """

    id: Optional[int] = Field(default=None, primary_key=True)
    project_id: int = Field(foreign_key="project.id")

    task_list: List["WorkflowTask"] = Relationship(
        sa_relationship_kwargs=dict(
            lazy="selectin",
            order_by="WorkflowTask.order",
            collection_class=ordering_list("order"),
            cascade="all, delete-orphan",
        ),
    )

    async def insert_task(
        self,
        task_id: int,
        *,
        args: Dict[str, Any] = None,
        order: Optional[int] = None,
        db: AsyncSession,
        commit: bool = True,
    ) -> WorkflowTask:
        if order is None:
            order = len(self.task_list)
        wf_task = WorkflowTask(task_id=task_id, args=args)
        db.add(wf_task)
        self.task_list.insert(order, wf_task)
        self.task_list.reorder()
        if commit:
            await db.commit()
        return wf_task
