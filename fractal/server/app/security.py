# Copyright 2022 (C)
# Original author: Jacopo Nespolo <jacopo.nespolo@exact-lab.it>
"""
Authentication and authorization services

The clients authenticate via LDAP or locally via OAUTH2
"""
from datetime import datetime
from datetime import timedelta
from datetime import timezone

from jose import jwt
from pydantic import BaseModel

from ..config import settings

# from fastapi.security import OAuth2PasswordBearer


# LDAPS
"""
DOC: http://www.pynut.com/?p=45
"""


class User(BaseModel):
    """
    Internal representation of a logged in user
    """

    sub: str


class Token(BaseModel):
    access_token: str
    token_type: str


class FailedAuthenticationException(Exception):
    pass


async def authenticate_user(username: str, password: str):
    """
    Authenticate user against LDAP

    Raise a FailedAuthenticationException if the user cannot authenticate
    """
    from ..utils import warn

    warn("authentication not implemented")
    # raise NotImplementedError
    return User(sub="1234")


def create_access_token(
    data: dict, expires_min: int = settings.JWT_EXPIRE_MINUTES
) -> str:
    """
    Create access token

    Paramters
    ---------
    data (dict): the payload
    expires_delta (timedelta): the lifetime of the token
    """
    payload = data.copy()
    expire = datetime.now(tz=timezone.utc) + timedelta(
        minutes=settings.JWT_EXPIRE_MINUTES
    )
    payload.update(dict(exp=expire))
    token = jwt.encode(
        payload, settings.JWT_SECRET_KEY, algorithm=settings.JWT_ALGORITHM
    )
    return token


# async def get_current_user(token: str = Depends(oauth2_scheme)):
#     pass
