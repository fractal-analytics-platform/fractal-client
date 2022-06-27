from pathlib import Path
from typing import Optional

from sqlmodel import Field
from sqlmodel import SQLModel


class Project(SQLModel, table=True):
    id: Optional[int] = Field(default=None, primary_key=True)
    name: str
    project_dir: Path

    # user_owner_id: UUID4 = Field(foreign_key="user_oauth.id")
    # user: Optional[UserOAuth] = Relationship(back_populates="oauth_accounts")
