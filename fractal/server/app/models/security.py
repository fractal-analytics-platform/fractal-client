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


def new_uuid() -> UUID4:
    # See https://github.com/tiangolo/sqlmodel/issues/25#issuecomment-982039809
    id_ = uuid.uuid4()
    while id_.hex[0] == "0":
        id_ = uuid.uuid4()
    return id_


class UserOAuth(SQLModelBaseUserDB, table=True):
    __tablename__ = "user_oauth"
    id: UUID4 = Field(
        default_factory=new_uuid,
        nullable=False,
        sa_column=Column(UUIDType(), primary_key=True),
    )
    oauth_accounts: List["OAuthAccount"] = Relationship(
        back_populates="user",
        sa_relationship_kwargs={"lazy": "selectin", "cascade": "all, delete"},
    )


class OAuthAccount(SQLModelBaseOAuthAccount, table=True):
    user_id: UUID4 = Field(foreign_key="user_oauth.id", nullable=False)
    user: Optional[UserOAuth] = Relationship(back_populates="oauth_accounts")


class UserRead(schemas.BaseUser[uuid.UUID]):
    pass


class UserUpdate(schemas.BaseUserUpdate):
    pass


class UserCreate(schemas.BaseUserCreate):
    pass
