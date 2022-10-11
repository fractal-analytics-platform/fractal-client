import uuid
from typing import List
from typing import Optional

from fastapi_users import schemas
from fastapi_users_db_sqlmodel import SQLModelBaseOAuthAccount
from fastapi_users_db_sqlmodel import SQLModelBaseUserDB
from pydantic import UUID4
from sqlalchemy_utils import UUIDType
from sqlmodel import Column
from sqlmodel import Field
from sqlmodel import Relationship


"""
Adapted from
    https://github.com/fastapi-users/fastapi-users-db-sqlmodel/
        blob/main/tests/conftest.py
"""


class UserOAuth(SQLModelBaseUserDB, table=True):
    __tablename__ = "user_oauth"
    id: UUID4 = Field(
        default_factory=uuid.uuid4,
        nullable=False,
        sa_column=Column(UUIDType(), primary_key=True),
    )
    slurm_user: Optional[str]
    oauth_accounts: List["OAuthAccount"] = Relationship(
        back_populates="user",
        sa_relationship_kwargs={"lazy": "selectin", "cascade": "all, delete"},
    )


class OAuthAccount(SQLModelBaseOAuthAccount, table=True):
    user_id: UUID4 = Field(foreign_key="user_oauth.id", nullable=False)
    user: Optional[UserOAuth] = Relationship(back_populates="oauth_accounts")


class UserRead(schemas.BaseUser[uuid.UUID]):
    slurm_user: str


class UserUpdate(schemas.BaseUserUpdate):
    pass


class UserCreate(schemas.BaseUserCreate):
    slurm_user: str
