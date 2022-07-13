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


class TaskBase(SQLModel):
    name: str = Field(sa_column_kwargs=dict(unique=True))
    resource_type: str
    module: Optional[str]
    input_type: str
    output_type: str
    default_parameters: Dict[str, Any] = Field(
        sa_column=Column(JSON), default={}
    )
    subtask_list: Optional[List["TaskBase"]] = Field(default=[])

    class Config:
        arbitrary_types_allowed = True


class TaskCreate(TaskBase):
    pass


class Subtask(SQLModel, table=True):  # type: ignore
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
    parameters: Dict[str, Any] = Field(sa_column=Column(JSON), default={})


class SubtaskRead(SQLModel):
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
    def callable(self):
        return self.module.partition(":")[2]

    @property
    def import_path(self):
        return self.module.partition(":")[0]


class TaskRead(TaskBase):
    id: int
    subtask_list: List[SubtaskRead]


SubtaskRead.update_forward_refs()
