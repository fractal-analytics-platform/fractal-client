# Copyright 2022 (C)
# Original author: Jacopo Nespolo <jacopo.nespolo@exact-lab.it>
"""
Authentication and authorization services

The clients authenticate via LDAP or locally via OAUTH2
"""
from datetime import datetime
from datetime import timedelta
from datetime import timezone
from typing import List

from fastapi import Depends
from fastapi.security import OAuth2PasswordBearer
from jose import jwt
from pydantic import BaseModel

from ..config import settings


# LDAPS
"""
DOC: http://www.pynut.com/?p=45
"""


oauth2_scheme = OAuth2PasswordBearer(
    tokenUrl="api/v1/token",
    scopes={
        "profile": "Can read information about self",
    },
)


class User(BaseModel):
    """
    Internal representation of a logged in user
    """

    sub: str
    name: str


class Token(BaseModel):
    access_token: str
    token_type: str


class TokenData(BaseModel):
    sub: str
    scopes: List[str] = []


class FailedAuthenticationException(Exception):
    pass


def mock_authenticate_user(username: str, password: str) -> User | None:
    from ..utils import warn

    warn("Using mock authentication")

    if username == "username" and password == "password":  # nosec
        return User(sub="1234", name="User Name")
    return None


def mock_get_user(sub: str) -> User:
    from ..utils import warn

    warn("Using mock get_user")

    if sub == "1234":
        return User(sub="1234", name="User Name")
    else:
        raise ValueError("User not found")


async def authenticate_user(username: str, password: str):
    """
    Authenticate user against LDAP

    Raise a FailedAuthenticationException if the user cannot authenticate
    """
    from ..utils import warn

    warn("Using mock authentication")

    user = mock_authenticate_user(username, password)
    if not user:
        raise FailedAuthenticationException
    return user


def access_token_encode(
    data: TokenData, expires_min: int = settings.JWT_EXPIRE_MINUTES
) -> str:
    """
    Create access token

    Paramters
    ---------
    data (dict): the payload
    expires_delta (timedelta): the lifetime of the token
    """
    payload = data.dict()
    expire = datetime.now(tz=timezone.utc) + timedelta(
        minutes=settings.JWT_EXPIRE_MINUTES
    )
    payload.update(dict(exp=expire))
    token = jwt.encode(
        payload, settings.JWT_SECRET_KEY, algorithm=settings.JWT_ALGORITHM
    )
    return token


def access_token_decode(token: str = Depends(oauth2_scheme)) -> User:
    payload = jwt.decode(token, settings.JWT_SECRET_KEY)
    sub = payload["sub"]
    return mock_get_user(sub)


async def token_authentication_login(
    username: str, password: str, scopes: List[str] | None = None
) -> Token:
    if not scopes:
        scopes = []
    user = await authenticate_user(username, password)
    access_token = access_token_encode(
        data=TokenData(sub=user.sub, scopes=scopes)
    )
    return Token(access_token=access_token, token_type="bearer")  # nosec


async def get_current_user(token: str = Depends(oauth2_scheme)):
    return access_token_decode(token)
