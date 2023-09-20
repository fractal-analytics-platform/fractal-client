from typing import Any
from typing import Optional

from pydantic import BaseModel
from pydantic import validator

from ._validators import valint
from ._validators import valstr
from .task import TaskExport
from .task import TaskImport
from .task import TaskRead


__all__ = (
    "WorkflowCreate",
    "WorkflowRead",
    "WorkflowUpdate",
    "WorkflowImport",
    "WorkflowExport",
    "WorkflowTaskCreate",
    "WorkflowTaskImport",
    "WorkflowTaskExport",
    "WorkflowTaskRead",
    "WorkflowTaskUpdate",
)


class _WorkflowTaskBase(BaseModel):

    meta: Optional[dict[str, Any]] = None
    args: Optional[dict[str, Any]] = None


class WorkflowTaskCreate(_WorkflowTaskBase):
    order: Optional[int]
    # Validators
    _order = validator("order", allow_reuse=True)(valint("order", min_val=0))


class WorkflowTaskRead(_WorkflowTaskBase):
    id: int
    order: Optional[int]
    workflow_id: int
    task_id: int
    task: TaskRead


class WorkflowTaskImport(_WorkflowTaskBase):
    task: TaskImport


class WorkflowTaskExport(_WorkflowTaskBase):
    task: TaskExport


class WorkflowTaskUpdate(_WorkflowTaskBase):
    # Validators
    @validator("meta")
    def check_no_parallelisation_level(cls, m):
        if "parallelization_level" in m:
            raise ValueError(
                "Overriding task parallelization level currently not allowed"
            )
        return m


class _WorkflowBase(BaseModel):
    name: str


class WorkflowRead(_WorkflowBase):
    id: int
    project_id: int
    task_list: list[WorkflowTaskRead]


class WorkflowCreate(_WorkflowBase):
    # Validators
    _name = validator("name", allow_reuse=True)(valstr("name"))


class WorkflowUpdate(_WorkflowBase):
    name: Optional[str]
    reordered_workflowtask_ids: Optional[list[int]]

    # Validators
    _name = validator("name", allow_reuse=True)(valstr("name"))

    @validator("reordered_workflowtask_ids")
    def check_positive_and_unique(cls, value):
        if any(i < 0 for i in value):
            raise ValueError("Negative `id` in `reordered_workflowtask_ids`")
        if len(value) != len(set(value)):
            raise ValueError("`reordered_workflowtask_ids` has repetitions")
        return value


class WorkflowImport(_WorkflowBase):
    task_list: list[WorkflowTaskImport]

    # Validators
    _name = validator("name", allow_reuse=True)(valstr("name"))


class WorkflowExport(_WorkflowBase):
    task_list: list[WorkflowTaskExport]
