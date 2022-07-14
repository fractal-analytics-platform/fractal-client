from enum import Enum
from typing import Any
from typing import Dict
from typing import List
from typing import Optional

from sqlalchemy import Column
from sqlalchemy.ext.asyncio import AsyncSession
from sqlalchemy.ext.orderinglist import ordering_list
from sqlalchemy.types import JSON
from sqlmodel import Field
from sqlmodel import Relationship
from sqlmodel import SQLModel


class ResourceTypeEnum(str, Enum):
    CORE_WORKFLOW = "core workflow"
    CORE_TASK = "core task"
    CORE_STEP = "core step"

    WORKFLOW = "workflow"
    TASK = "task"
    STEP = "step"


class TaskBase(SQLModel):
    name: str = Field(sa_column_kwargs=dict(unique=True))
    resource_type: ResourceTypeEnum
    module: Optional[str]
    input_type: str
    output_type: str
    default_args: Dict[str, Any] = Field(sa_column=Column(JSON), default={})
    subtask_list: Optional[List["TaskBase"]] = Field(default=[])

    class Config:
        arbitrary_types_allowed = True


class TaskCreate(TaskBase):
    pass


class SubtaskBase(SQLModel):
    pass


class Subtask(SubtaskBase, table=True):  # type: ignore
    class Config:
        arbitrary_types_allowed = True
        fields = {"parent": {"exclude": True}}

    parent_task_id: Optional[int] = Field(
        default=None, foreign_key="task.id", primary_key=True
    )
    subtask_id: Optional[int] = Field(
        default=None, foreign_key="task.id", primary_key=True
    )
    parent: "Task" = Relationship(
        sa_relationship_kwargs=dict(
            primaryjoin="Subtask.parent_task_id==Task.id"
        )
    )
    subtask: "Task" = Relationship(
        sa_relationship_kwargs=dict(primaryjoin="Subtask.subtask_id==Task.id")
    )
    order: Optional[int]
    args: Dict[str, Any] = Field(sa_column=Column(JSON), default={})

    @property
    def _is_atomic(self):
        """
        A subtask is atomic iff the child task is atomic
        """
        return self.subtask._is_atomic

    @property
    def _merged_args(self):
        out = self.subtask.default_args.copy()
        out.update(self.args)
        return out

    @property
    def import_path(self):
        return self.subtask.import_path

    @property
    def callable(self):
        return self.subtask.callable


class SubtaskRead(SubtaskBase):
    subtask: "TaskRead"


class Task(TaskBase, table=True):  # type: ignore
    id: Optional[int] = Field(default=None, primary_key=True)
    subtask_list: List["Subtask"] = Relationship(
        sa_relationship_kwargs=dict(
            primaryjoin="Task.id==Subtask.parent_task_id",
            lazy="selectin",
            order_by="Subtask.order",
            collection_class=ordering_list("order"),
        ),
    )

    @property
    def _is_atomic(self) -> bool:
        """
        A task is atomic iff it does not contain subtasks
        """
        return bool(self.subtask_list)

    @property
    def callable(self):
        return self.module.partition(":")[2]

    @property
    def import_path(self):
        return self.module.partition(":")[0]

    async def add_subtask(
        self,
        db: AsyncSession,
        subtask: "Task",
        order=None,
        args=None,
        commit_and_refresh=True,
    ):
        if not args:
            args = dict()
        st = Subtask(parent=self, subtask=subtask, args=args)

        if order is None:
            self.subtask_list.append(st)
        else:
            self.subtask_list.insert(order, st)
        db.add_all([self, st])
        if commit_and_refresh:
            await db.commit()
            await db.refresh(self)


class TaskRead(TaskBase):
    id: int
    subtask_list: List[SubtaskRead]


SubtaskRead.update_forward_refs()
