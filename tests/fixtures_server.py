import asyncio
from dataclasses import dataclass
from typing import Any
from typing import AsyncGenerator
from typing import List
from typing import Optional

import pytest
from asgi_lifespan import LifespanManager
from fastapi import FastAPI
from httpx import AsyncClient
from sqlalchemy.ext.asyncio import AsyncEngine
from sqlalchemy.ext.asyncio import AsyncSession
from sqlalchemy.ext.asyncio import create_async_engine
from sqlalchemy.orm import sessionmaker

from fractal.server import start_application


@pytest.fixture(scope="session")
def event_loop():
    _event_loop = asyncio.new_event_loop()
    yield _event_loop


@pytest.fixture(scope="session")
def patch_settings():
    from os import environ

    environ["JWT_SECRET_KEY"] = "secret_key"
    environ["DEPLOYMENT_TYPE"] = "development"
    environ["DATA_DIR_ROOT"] = "/tmp/"

    environ["DB_ENGINE"] = "sqlite"
    environ["SQLITE_PATH"] = ""  # in memory

    from fractal.server.config import settings

    settings.LDAP_SERVER = "testserver"
    return settings


@pytest.fixture(scope="session")
async def db_engine(patch_settings) -> AsyncGenerator[AsyncEngine, None]:
    from fractal.server.config import settings

    engine = create_async_engine(
        settings.DATABASE_URL, echo=False, future=True
    )
    yield engine


@pytest.fixture
async def db(db_engine, app) -> AsyncGenerator[AsyncSession, None]:
    from sqlmodel import SQLModel

    async with db_engine.begin() as conn:
        await conn.run_sync(SQLModel.metadata.create_all)
        async_session_maker = sessionmaker(
            db_engine, class_=AsyncSession, expire_on_commit=False
        )
        async with async_session_maker() as session:
            from fractal.server.app.db import get_db

            app.dependency_overrides[get_db] = lambda: session
            yield session
        await conn.run_sync(SQLModel.metadata.drop_all)


@pytest.fixture
async def app(patch_settings) -> AsyncGenerator[FastAPI, Any]:
    yield start_application()


@pytest.fixture
async def client(app: FastAPI) -> AsyncGenerator[AsyncClient, Any]:
    async with AsyncClient(
        app=app, base_url="http://test"
    ) as client, LifespanManager(app):
        yield client


@pytest.fixture
async def MockCurrentUser(app):
    from fractal.server.app.security import current_active_user
    from fractal.server.app.security import User

    @dataclass
    class _MockCurrentUser:
        """
        Context managed user override
        """

        sub: str
        scopes: Optional[List[str]]

        def current_active_user_override(self):
            def __current_active_user_override():
                return User(sub=self.sub, name=self.sub)

            return __current_active_user_override

        def __enter__(self):
            self.previous_user = app.dependency_overrides.get(
                current_active_user, None
            )
            app.dependency_overrides[
                current_active_user
            ] = self.current_active_user_override()

        def __exit__(self, *args, **kwargs):
            if self.previous_user:
                app.dependency_overrides[
                    current_active_user
                ] = self.previous_user()

    return _MockCurrentUser
