from typing import Any
from typing import Dict
from typing import List
from typing import Optional

from pydantic import validator
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

    @validator("args")
    def validate_args(cls, value):
        forbidden_args_keys = {
            "input_paths",
            "output_path",
            "metadata",
            "logger",
        }
        args_keys = set(value.keys())
        intersect_keys = forbidden_args_keys.intersection(args_keys)
        if intersect_keys:
            raise ValueError(
                "`args` contains the following forbidden keys: "
                f"{intersect_keys}"
            )
        return value

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

    @property
    def is_parallel(self) -> bool:
        return self.task.is_parallel

    @property
    def parallelization_level(self) -> str:
        return self.task.parallelization_level

    def assemble_args(self, extra: Dict[str, Any] = None):
        """
        Merge of `extra` arguments and `self.arguments`.

        Return
        ------
        full_arsgs (Dict):
            A dictionary consisting of the merge of `extra` and
            self.arguments.
        """
        full_args = {}
        if extra:
            full_args.update(extra)
        full_args.update(self.arguments)
        return full_args


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
