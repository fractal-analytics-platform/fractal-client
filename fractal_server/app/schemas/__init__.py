from datetime import datetime
from enum import Enum
from pathlib import Path
from typing import Any
from typing import Dict
from typing import List
from typing import Optional

from pydantic import root_validator
from pydantic import validator
from sqlmodel import Field
from sqlmodel import SQLModel

from ...utils import slugify


__all__ = (
    "ApplyWorkflowBase",
    "ApplyWorkflowCreate",
    "ApplyWorkflowRead",
    "ProjectBase",
    "ProjectCreate",
    "ProjectRead",
    "DatasetBase",
    "DatasetUpdate",
    "DatasetCreate",
    "DatasetRead",
    "ResourceBase",
    "ResourceCreate",
    "ResourceRead",
    "ResourceUpdate",
    "TaskBase",
    "TaskCreate",
    "TaskUpdate",
    "TaskRead",
)


class ApplyWorkflowBase(SQLModel):
    project_id: int
    input_dataset_id: int
    output_dataset_id: Optional[int]
    workflow_id: Optional[int]
    overwrite_input: bool = False


class ApplyWorkflowCreate(ApplyWorkflowBase):
    pass


class ApplyWorkflowRead(ApplyWorkflowBase):
    id: int
    start_timestamp: datetime


# PROJECT


class ProjectBase(SQLModel):
    name: str
    project_dir: str
    read_only: bool = False


class ProjectCreate(ProjectBase):
    slug: Optional[str] = Field()
    default_dataset_name: Optional[str] = "default"

    @root_validator(pre=True)
    def compute_slug_if_not_provided(cls, values):
        """
        Creates the project's slug

        Or double checks the user provided one if any.
        """
        slug = values.get("slug", None)
        if not slug:
            slug = values["name"]
            values["slug"] = slugify(slug)
        return values

    @validator("default_dataset_name")
    def not_null(cls, value):
        if not value:
            value = "default"
        return value


class ProjectRead(ProjectBase):
    id: int
    slug: str
    dataset_list: List["DatasetRead"] = []


# DATASET


class DatasetBase(SQLModel):
    name: str
    project_id: Optional[int]
    type: Optional[str]
    meta: Dict[str, Any] = {}
    read_only: Optional[bool] = False


class DatasetUpdate(DatasetBase):
    name: Optional[str]  # type:ignore
    meta: Optional[Dict[str, Any]] = None  # type:ignore


class DatasetCreate(DatasetBase):
    pass


class DatasetRead(DatasetBase):
    id: int
    resource_list: List["ResourceRead"]


# RESOURCE


class ResourceBase(SQLModel):
    path: str
    glob_pattern: Optional[str] = ""

    @property
    def glob_path(self) -> Path:
        return Path(self.path) / self.glob_pattern


class ResourceCreate(ResourceBase):
    pass


class ResourceUpdate(ResourceBase):
    pass


class ResourceRead(ResourceBase):
    id: int
    dataset_id: int


ProjectRead.update_forward_refs()


# TASK


class ResourceTypeEnum(str, Enum):
    CORE_WORKFLOW = "core workflow"
    CORE_TASK = "core task"
    CORE_STEP = "core step"

    WORKFLOW = "workflow"
    TASK = "task"
    STEP = "step"


class TaskBase(SQLModel):
    name: str
    resource_type: ResourceTypeEnum
    module: Optional[str]
    input_type: str
    output_type: str
    default_args: Dict[str, Any] = Field(default={})
    project_id: Optional[int]

    class Config:
        arbitrary_types_allowed = True


class TaskUpdate(TaskBase):
    name: Optional[str]  # type:ignore
    input_type: Optional[str]  # type:ignore
    output_type: Optional[str]  # type:ignore
    default_args: Optional[Dict[str, Any]] = None  # type:ignore


class TaskCreate(TaskBase):
    pass


class TaskRead(TaskBase):
    id: int


DatasetRead.update_forward_refs()
