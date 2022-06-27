from typing import List
from typing import Optional

from fastapi_users_db_sqlmodel import SQLModelBaseOAuthAccount
from fastapi_users_db_sqlmodel import SQLModelBaseUserDB
from pydantic import BaseModel
from pydantic import UUID4
from sqlmodel import Field
from sqlmodel import Relationship


"""
Adapted from
    https://github.com/fastapi-users/fastapi-users-db-sqlmodel/
        blob/main/tests/conftest.py
"""


class UserOAuth(SQLModelBaseUserDB, table=True):
    __tablename__ = "user_oauth"
    oauth_accounts: List["OAuthAccount"] = Relationship(
        back_populates="user",
        sa_relationship_kwargs={"lazy": "selectin", "cascade": "all, delete"},
    )


class OAuthAccount(SQLModelBaseOAuthAccount, table=True):
    user_id: UUID4 = Field(foreign_key="user_oauth.id")
    user: Optional[UserOAuth] = Relationship(back_populates="oauth_accounts")


class UserRead(BaseModel):
    pass


class UserUpdate(BaseModel):
    pass
