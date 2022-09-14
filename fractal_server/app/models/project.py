from pathlib import Path
from typing import Any
from typing import Dict
from typing import List
from typing import Optional

from fractal.common.models import DatasetBase
from fractal.common.models import ProjectBase
from fractal.common.models import ResourceBase
from pydantic import UUID4
from sqlalchemy import Column
from sqlalchemy.types import JSON
from sqlmodel import Field
from sqlmodel import Relationship

from .security import UserOAuth


class Dataset(DatasetBase, table=True):  # type: ignore
    id: Optional[int] = Field(default=None, primary_key=True)
    project_id: int = Field(foreign_key="project.id")
    resource_list: List["Resource"] = Relationship(
        sa_relationship_kwargs={"lazy": "selectin"}
    )
    meta: Dict[str, Any] = Field(sa_column=Column(JSON), default={})

    class Config:
        arbitrary_types_allowed = True

    @property
    def paths(self) -> List[Path]:
        return [r.glob_path for r in self.resource_list]


class Project(ProjectBase, table=True):  # type: ignore
    id: Optional[int] = Field(default=None, primary_key=True)
    slug: Optional[str] = Field()

    user_owner_id: Optional[UUID4] = Field(foreign_key="user_oauth.id")
    user: Optional[UserOAuth] = Relationship()

    dataset_list: List[Dataset] = Relationship(
        sa_relationship_kwargs={
            "lazy": "selectin",
            "cascade": "all, delete-orphan",
        }
    )


class Resource(ResourceBase, table=True):  # type: ignore
    id: Optional[int] = Field(default=None, primary_key=True)
    dataset_id: int = Field(foreign_key="dataset.id")
