import uuid
from typing import AsyncGenerator

from fastapi import APIRouter
from fastapi import Depends
from fastapi_users import BaseUserManager
from fastapi_users import FastAPIUsers
from fastapi_users import UUIDIDMixin
from fastapi_users.authentication import AuthenticationBackend
from fastapi_users.authentication import BearerTransport
from fastapi_users.authentication import CookieTransport
from fastapi_users.authentication import JWTStrategy
from fastapi_users_db_sqlmodel import SQLModelUserDatabaseAsync
from httpx_oauth.clients.github import GitHubOAuth2
from sqlalchemy.ext.asyncio import AsyncSession

from ...config import settings
from ..db import get_db
from ..models.security import OAuthAccount
from ..models.security import UserCreate
from ..models.security import UserOAuth as User
from ..models.security import UserRead
from ..models.security import UserUpdate


github_client = GitHubOAuth2(
    settings.OAUTH_GITHUB_CLIENT_ID, settings.OAUTH_GITHUB_CLIENT_SECRET
)


async def get_user_db(
    session: AsyncSession = Depends(get_db),
) -> SQLModelUserDatabaseAsync:
    yield SQLModelUserDatabaseAsync(session, User, OAuthAccount)


class UserManager(UUIDIDMixin, BaseUserManager[User, uuid.UUID]):
    pass


async def get_user_manager(
    user_db: SQLModelUserDatabaseAsync = Depends(get_user_db),
) -> AsyncGenerator[UserManager, None]:
    yield UserManager(user_db)


bearer_transport = BearerTransport(tokenUrl="/auth/token/login")
cookie_transport = CookieTransport()


def get_jwt_strategy() -> JWTStrategy:
    return JWTStrategy(
        secret=settings.JWT_SECRET_KEY,
        lifetime_seconds=settings.JWT_EXPIRE_SECONDS,
    )


token_backend = AuthenticationBackend(
    name="bearer-jwt",
    transport=bearer_transport,
    get_strategy=get_jwt_strategy,
)
cookie_backend = AuthenticationBackend(
    name="cookie-jwt",
    transport=cookie_transport,
    get_strategy=get_jwt_strategy,
)


fastapi_users = FastAPIUsers[User, uuid.UUID](
    get_user_manager,
    [token_backend, cookie_backend],
)


current_active_user = fastapi_users.current_user(active=True)


# AUTH ROUTES

auth_router = APIRouter()

auth_router.include_router(
    fastapi_users.get_auth_router(token_backend),
    prefix="/token",
)
auth_router.include_router(
    fastapi_users.get_auth_router(cookie_backend),
)
auth_router.include_router(
    fastapi_users.get_register_router(UserRead, UserCreate),
)
auth_router.include_router(
    fastapi_users.get_reset_password_router(),
)
auth_router.include_router(
    fastapi_users.get_verify_router(UserRead),
)
auth_router.include_router(
    fastapi_users.get_users_router(UserRead, UserUpdate),
    prefix="/users",
)

# GitHub OAuth
auth_router.include_router(
    fastapi_users.get_oauth_router(
        github_client,
        cookie_backend,
        settings.JWT_SECRET_KEY,
        # WARNING:
        # associate_by_email=True exposes to security risks if the OAuth
        # provider does not verify emails.
        associate_by_email=True,
    ),
    prefix="/github",
)
auth_router.include_router(
    fastapi_users.get_oauth_associate_router(
        github_client, UserRead, settings.JWT_SECRET_KEY
    ),
    prefix="/github/associate",
)
