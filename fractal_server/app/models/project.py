from pathlib import Path
from typing import Any
from typing import Dict
from typing import List
from typing import Optional

from pydantic import UUID4
from sqlalchemy import Column
from sqlalchemy.types import JSON
from sqlmodel import Field
from sqlmodel import Relationship
from sqlmodel import SQLModel

from ..schemas.project import _DatasetBase
from ..schemas.project import _ProjectBase
from ..schemas.project import _ResourceBase
from .security import UserOAuth as User


class LinkUserProject(SQLModel, table=True):
    """
    Crossing table between User and Project
    """

    project_id: int = Field(foreign_key="project.id", primary_key=True)
    user_id: UUID4 = Field(foreign_key="user_oauth.id", primary_key=True)


class Dataset(_DatasetBase, table=True):
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


class Project(_ProjectBase, table=True):
    id: Optional[int] = Field(default=None, primary_key=True)
    slug: Optional[str] = Field()

    user_member_list: List[User] = Relationship(
        link_model=LinkUserProject,
        sa_relationship_kwargs={
            "lazy": "selectin",
        },
    )

    dataset_list: List[Dataset] = Relationship(
        sa_relationship_kwargs={
            "lazy": "selectin",
            "cascade": "all, delete-orphan",
        }
    )


class Resource(_ResourceBase, table=True):
    id: Optional[int] = Field(default=None, primary_key=True)
    dataset_id: int = Field(foreign_key="dataset.id")
