from typing import Optional

from pydantic import root_validator
from pydantic import UUID4
from sqlmodel import Field
from sqlmodel import Relationship
from sqlmodel import SQLModel

from .security import UserOAuth


def slugify(value: str):
    return value.lower().replace(" ", "_")


class ProjectBase(SQLModel):
    name: str
    project_dir: str


class ProjectCreate(ProjectBase):
    slug: Optional[str] = Field(sa_column_kwargs={"unique": True})

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


class ProjectRead(ProjectBase):
    id: int
    slug: str
