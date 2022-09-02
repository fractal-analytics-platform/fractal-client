from enum import Enum
from pathlib import Path
from typing import Any
from typing import Dict
from typing import List
from typing import Optional

from pydantic import BaseModel
from pydantic import root_validator
from sqlmodel import Field
from sqlmodel import SQLModel

from ..utils import slugify


class ApplyWorkflow(BaseModel):
    project_id: int
    input_dataset_id: int
    output_dataset_id: Optional[int]
    workflow_id: Optional[int]
    overwrite_input: bool = False


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
    name: str = Field(sa_column_kwargs=dict(unique=True))
    resource_type: ResourceTypeEnum
    module: Optional[str]
    input_type: str
    output_type: str
    default_args: Dict[str, Any] = Field(default={})

    class Config:
        arbitrary_types_allowed = True


class TaskUpdate(TaskBase):
    name: Optional[str]  # type:ignore
    resource_type: Optional[ResourceTypeEnum]  # type:ignore
    input_type: Optional[str]  # type:ignore
    output_type: Optional[str]  # type:ignore
    default_args: Optional[Dict[str, Any]] = None  # type:ignore
    subtask_list: Optional[List["TaskBase"]] = Field(default=[])


class TaskCreate(TaskBase):
    pass


class SubtaskBase(SQLModel):
    parent_task_id: Optional[int] = None
    subtask_id: Optional[int] = None
    order: Optional[int] = None
    args: Dict[str, Any] = Field(default={})


class SubtaskCreate(SubtaskBase):
    subtask_id: int


class SubtaskRead(SubtaskBase):
    parent_task_id: int
    subtask_id: int
    subtask: "TaskRead"


class TaskRead(TaskBase):
    id: int
    subtask_list: List[SubtaskRead]


SubtaskRead.update_forward_refs()
DatasetRead.update_forward_refs()
