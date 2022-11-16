import uuid
from typing import AsyncGenerator
from typing import Optional
from typing import Union

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
from sqlalchemy.ext.asyncio import AsyncSession

from ...config import get_settings
from ...syringe import Inject
from ..db import get_db
from ..models.security import OAuthAccount
from ..models.security import UserCreate
from ..models.security import UserOAuth as User
from ..models.security import UserRead
from ..models.security import UserUpdate


async def get_user_db(
    session: AsyncSession = Depends(get_db),
) -> AsyncGenerator[SQLModelUserDatabaseAsync, None]:
    yield SQLModelUserDatabaseAsync(session, User, OAuthAccount)


class UserManager(UUIDIDMixin, BaseUserManager[User, uuid.UUID]):
    pass


async def get_user_manager(
    user_db: SQLModelUserDatabaseAsync = Depends(get_user_db),
) -> AsyncGenerator[UserManager, None]:
    yield UserManager(user_db)


bearer_transport = BearerTransport(tokenUrl="/auth/token/login")
cookie_transport = CookieTransport(cookie_samesite="none")


def get_jwt_strategy() -> JWTStrategy:
    settings = Inject(get_settings)
    return JWTStrategy(
        secret=settings.JWT_SECRET_KEY,
        lifetime_seconds=settings.JWT_EXPIRE_SECONDS,
    )


def get_jwt_cookie_strategy() -> JWTStrategy:
    settings = Inject(get_settings)
    return JWTStrategy(
        secret=settings.JWT_SECRET_KEY,
        lifetime_seconds=settings.COOKIE_EXPIRE_SECONDS,
    )


token_backend = AuthenticationBackend(
    name="bearer-jwt",
    transport=bearer_transport,
    get_strategy=get_jwt_strategy,
)
cookie_backend = AuthenticationBackend(
    name="cookie-jwt",
    transport=cookie_transport,
    get_strategy=get_jwt_cookie_strategy,
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


# OAUTH CLIENTS
settings = Inject(get_settings)
for client in settings.OAUTH_CLIENTS:
    # INIT CLIENTS
    client_name = client.CLIENT_NAME.lower()
    _client: Optional[Union["GitHubOAuth2", "OAuth2"]] = None
    if client_name == "github":
        from httpx_oauth.clients.github import GitHubOAuth2

        _client = GitHubOAuth2(client.CLIENT_ID, client.CLIENT_SECRET)
    else:  # GENERIC CLIENT
        from httpx_oauth.oauth2 import OAuth2

        if (
            not client.CLIENT_SECRET
            or not client.AUTHORIZE_ENDPOINT
            or not client.ACCESS_TOKEN_ENDPOINT
        ):
            raise ValueError(
                "Must specify CLIENT_SECRET, AUTHORIZE_ENDPOINT and "
                "ACCESS_TOKEN_ENDPOINT to define custom OAuth2 client."
            )
        _client = OAuth2(
            client.CLIENT_ID,
            client.CLIENT_SECRET,
            client.AUTHORIZE_ENDPOINT,
            client.ACCESS_TOKEN_ENDPOINT,
            refresh_token_endpoint=client.REFRESH_TOKEN_ENDPOINT,
            revoke_token_endpoint=client.REVOKE_TOKEN_ENDPOINT,
        )

    # ADD ROUTES
    # GitHub OAuth
    auth_router.include_router(
        fastapi_users.get_oauth_router(
            _client,
            cookie_backend,
            settings.JWT_SECRET_KEY,
            # WARNING:
            # associate_by_email=True exposes to security risks if the OAuth
            # provider does not verify emails.
            associate_by_email=True,
        ),
        prefix=f"/{client_name}",
    )
    auth_router.include_router(
        fastapi_users.get_oauth_associate_router(
            _client, UserRead, settings.JWT_SECRET_KEY
        ),
        prefix=f"/{client_name}/associate",
    )
