from datetime import datetime
from enum import Enum
from pathlib import Path
from typing import Optional

from sqlalchemy import Column
from sqlalchemy.types import DateTime
from sqlmodel import Field
from sqlmodel import Relationship

from ...config import get_settings
from ...syringe import Inject
from ...utils import get_timestamp
from ..schemas import ApplyWorkflowBase
from .project import Dataset
from .project import Project
from .workflow import Workflow


class StatusType(str, Enum):
    SUBMITTED = "submitted"
    PENDING = "pending"
    DONE = "done"


class ApplyWorkflow(ApplyWorkflowBase, table=True):  # type: ignore
    id: Optional[int] = Field(default=None, primary_key=True)
    project_id: int = Field(foreign_key="project.id")
    input_dataset_id: int = Field(foreign_key="dataset.id")
    output_dataset_id: int = Field(foreign_key="dataset.id")
    workflow_id: int = Field(foreign_key="workflow.id")

    project: Project = Relationship()
    input_dataset: Dataset = Relationship(
        sa_relationship_kwargs=dict(
            lazy="selectin",
            primaryjoin="ApplyWorkflow.input_dataset_id==Dataset.id",
        )
    )
    output_dataset: Dataset = Relationship(
        sa_relationship_kwargs=dict(
            lazy="selectin",
            primaryjoin="ApplyWorkflow.output_dataset_id==Dataset.id",
        )
    )
    workflow: Workflow = Relationship()

    start_timestamp: datetime = Field(
        default_factory=get_timestamp,
        nullable=False,
        sa_column=Column(DateTime(timezone=True)),
    )
    status: StatusType = StatusType.SUBMITTED

    @property
    def job_root_path(self) -> Path:
        settings = Inject(get_settings)
        return settings.RUNNER_ROOT_DIR / f"job_{self.id:06d}"

    @property
    def log_path(self) -> Path:
        return self.job_root_path / "job.log"

    def make_job_dir(self):
        _path = self.job_root_path
        if not _path.is_dir():
            _path.mkdir(exists_ok=True, parents=True)
