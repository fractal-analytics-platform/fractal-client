from dataclasses import dataclass
from typing import Any
from typing import AsyncGenerator
from typing import List

import pytest
from asgi_lifespan import LifespanManager
from fastapi import FastAPI
from httpx import AsyncClient

from fractal.server import start_application


@pytest.fixture(scope="session")
def patch_settings():
    from os import environ

    environ["JWT_SECRET_KEY"] = "secret_key"
    environ["DEPLOYMENT_TYPE"] = "development"
    environ["DATA_DIR_ROOT"] = "/tmp/"

    from fractal.server.config import settings

    settings.LDAP_SERVER = "testserver"
    return settings


@pytest.fixture
async def app(patch_settings) -> AsyncGenerator[FastAPI, Any]:
    yield start_application()


@pytest.fixture
async def client(app: FastAPI) -> AsyncGenerator[AsyncClient, Any]:
    async with AsyncClient(
        app=app, base_url="http://test/api"
    ) as client, LifespanManager(app):
        yield client


@pytest.fixture
async def mock_ldap(patch_settings):
    from ldap3 import Server

    server = Server(patch_settings.LDAP_SERVER)
    yield server


@pytest.fixture
async def MockCurrentUser(patch_settings):
    from fractal.server.app.security import get_current_user
    from fractal.server.app.security import User

    @dataclass
    class _MockCurrentUser:
        """
        Context managed user override
        """

        app: FastAPI
        sub: str
        scopes: List[str] | None

        def get_user_override(self):
            def __get_user_override():
                return User(sub=self.sub, name=self.sub)

        def __enter__(self):
            self.previous_user = self.app.dependency_overrides.get(
                get_current_user, None
            )
            self.app.dependency_overrides[
                get_current_user
            ] = self.get_user_override()

        def __exit__(self, *args, **kwargs):
            if self.previous_user:
                self.app.dependency_overrides[
                    get_current_user
                ] = self.previous_user()

    return _MockCurrentUser
