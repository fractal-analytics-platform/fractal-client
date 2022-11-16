from datetime import datetime
from typing import Optional

from sqlmodel import SQLModel

from .manifest import *  # noqa: F403
from .project import *  # noqa: F403
from .state import *  # noqa: F403
from .task import *  # noqa: F403
from .workflow import *  # noqa: F403


__all__ = (
    (
        "ApplyWorkflowBase",
        "ApplyWorkflowCreate",
        "ApplyWorkflowRead",
    )
    + project.__all__  # noqa: F405
    + task.__all__  # noqa: F405
    + workflow.__all__  # noqa: F405
    + manifest.__all__  # noqa: F405
    + state.__all__  # noqa: F405
)


class ApplyWorkflowBase(SQLModel):
    project_id: int
    input_dataset_id: int
    output_dataset_id: Optional[int]
    workflow_id: Optional[int]
    overwrite_input: bool = False
    worker_init: Optional[str]


class ApplyWorkflowCreate(ApplyWorkflowBase):
    pass


class ApplyWorkflowRead(ApplyWorkflowBase):
    id: int
    start_timestamp: datetime
