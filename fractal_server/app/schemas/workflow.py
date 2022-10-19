from typing import List

from sqlmodel import SQLModel

from .task import TaskRead


__all__ = ("WorkflowRead", "WorkflowCreate")


class _WorkflowBase(SQLModel):
    name: str
    project_id: int


class WorkflowRead(_WorkflowBase):
    id: int
    task_list: List[TaskRead]


class WorkflowCreate(_WorkflowBase):
    pass
