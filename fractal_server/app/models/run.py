from datetime import datetime
from enum import Enum
from typing import Optional

from fractal.common.models import ApplyWorkflowBase
from sqlmodel import Field

from .models_utils import get_timestamp


class StatusType(str, Enum):
    SUBMITTED = "submitted"
    PENDING = "pending"
    DONE = "done"


class ApplyWorkflow(ApplyWorkflowBase, table=True):  # type: ignore
    id: Optional[int] = Field(default=None, primary_key=True)
    project_id: int = Field(foreign_key="project.id")
    input_dataset_id: int = Field(foreign_key="dataset.id")
    output_dataset_id: int = Field(foreign_key="dataset.id")
    workflow_id: int = Field(foreign_key="task.id")

    start_timestamp: datetime = Field(
        default_factory=get_timestamp, nullable=False
    )
    status: StatusType = StatusType.SUBMITTED
