from pathlib import Path
from typing import Any
from typing import Dict
from typing import List
from typing import Optional

from pydantic import validator
from sqlmodel import SQLModel

__all__ = (
    "ProjectCreate",
    "ProjectRead",
    "DatasetUpdate",
    "DatasetCreate",
    "DatasetRead",
    "ResourceCreate",
    "ResourceRead",
    "ResourceUpdate",
)


class _ProjectBase(SQLModel):
    name: str
    project_dir: str
    read_only: bool = False


class ProjectCreate(_ProjectBase):
    default_dataset_name: Optional[str] = "default"

    @validator("default_dataset_name")
    def not_null(cls, value):
        if not value:
            value = "default"
        return value


class ProjectRead(_ProjectBase):
    id: int
    slug: str
    dataset_list: List["DatasetRead"] = []


# DATASET


class _DatasetBase(SQLModel):
    name: str
    project_id: Optional[int]
    type: Optional[str]
    meta: Dict[str, Any] = {}
    read_only: Optional[bool] = False


class DatasetUpdate(_DatasetBase):
    name: Optional[str]  # type:ignore
    meta: Optional[Dict[str, Any]] = None  # type:ignore


class DatasetCreate(_DatasetBase):
    pass


class DatasetRead(_DatasetBase):
    id: int
    resource_list: List["ResourceRead"]


# RESOURCE


class _ResourceBase(SQLModel):
    path: str
    glob_pattern: Optional[str] = ""

    @property
    def glob_path(self) -> Path:
        return Path(self.path) / self.glob_pattern


class ResourceCreate(_ResourceBase):
    pass


class ResourceUpdate(_ResourceBase):
    pass


class ResourceRead(_ResourceBase):
    id: int
    dataset_id: int


ProjectRead.update_forward_refs()
DatasetRead.update_forward_refs()
