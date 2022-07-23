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
    slug: Optional[str] = Field(sa_column_kwargs={"unique": True})
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


class DatasetCreate(DatasetBase):
    pass


class DatasetRead(DatasetBase):
    id: int


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
