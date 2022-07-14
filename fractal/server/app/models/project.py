from pathlib import Path
from typing import Any
from typing import Dict
from typing import List
from typing import Optional

from pydantic import root_validator
from pydantic import UUID4
from sqlalchemy import Column
from sqlalchemy.types import JSON
from sqlmodel import Field
from sqlmodel import Relationship
from sqlmodel import SQLModel

from .security import UserOAuth


def slugify(value: str):
    return value.lower().replace(" ", "_")


class DatasetBase(SQLModel):
    name: str
    type: Optional[str]
    meta: Dict[str, Any] = {}


class DatasetCreate(DatasetBase):
    pass


class DatasetRead(DatasetBase):
    id: int


class ProjectBase(SQLModel):
    name: str
    project_dir: str
    read_only: bool = False


class ResourceBase(SQLModel):
    path: str
    glob_pattern: str

    @property
    def glob_path(self):
        return Path(self.path) / self.glob_pattern


class Dataset(DatasetBase, table=True):  # type: ignore
    id: Optional[int] = Field(default=None, primary_key=True)
    project_id: int = Field(foreign_key="project.id")
    resource_list: List["Resource"] = Relationship(
        sa_relationship_kwargs={"lazy": "selectin"}
    )
    meta: Dict[str, Any] = Field(sa_column=Column(JSON), default={})

    class Config:
        arbitrary_types_allowed = True


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


class Project(ProjectBase, table=True):  # type: ignore
    id: Optional[int] = Field(default=None, primary_key=True)
    slug: str

    user_owner_id: Optional[UUID4] = Field(foreign_key="user_oauth.id")
    user: Optional[UserOAuth] = Relationship()

    dataset_list: List[Dataset] = Relationship(
        sa_relationship_kwargs={"lazy": "selectin"}
    )


class ResourceCreate(ResourceBase):
    pass


class Resource(ResourceBase, table=True):  # type: ignore
    id: Optional[int] = Field(default=None, primary_key=True)
    dataset_id: int = Field(foreign_key="dataset.id")


class ProjectRead(ProjectBase):
    id: int
    slug: str
    dataset_list: List[DatasetRead] = []


class ResourceRead(ResourceBase):
    id: int
    dataset_id: int
