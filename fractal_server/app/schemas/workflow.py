from typing import Any
from typing import Dict
from typing import List
from typing import Optional

from sqlmodel import SQLModel

from .task import TaskRead


__all__ = (
    "WorkflowCreate",
    "WorkflowRead",
    "WorkflowUpdate",
    "WorkflowTaskCreate",
    "WorkflowTaskRead",
    "WorkflowTaskUpdate",
)


class _WorkflowTaskBase(SQLModel):
    workflow_id: Optional[int]
    task_id: Optional[int]
    order: Optional[int]
    args: Optional[Dict[str, Any]] = None


class WorkflowTaskCreate(_WorkflowTaskBase):
    task_id: int


class WorkflowTaskRead(_WorkflowTaskBase):
    id: int
    workflow_id: int
    task: TaskRead


class WorkflowTaskUpdate(_WorkflowTaskBase):
    args: Optional[Dict[str, Any]]  # type: ignore


class _WorkflowBase(SQLModel):
    name: str
    project_id: int


class WorkflowRead(_WorkflowBase):
    id: int
    task_list: List[WorkflowTaskRead]


class WorkflowCreate(_WorkflowBase):
    pass


class WorkflowUpdate(_WorkflowBase):
    name: Optional[str]  # type: ignore
    project_id: Optional[int]  # type: ignore
