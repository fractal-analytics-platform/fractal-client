from typing import Any
from typing import AsyncGenerator

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
