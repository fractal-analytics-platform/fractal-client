from typing import Any
from typing import AsyncGenerator

import pytest
from asgi_lifespan import LifespanManager
from fastapi import FastAPI
from httpx import AsyncClient

from fractal.server import start_application


@pytest.fixture
async def app() -> AsyncGenerator[FastAPI, Any]:
    yield start_application()


@pytest.fixture
async def client(app: FastAPI) -> AsyncGenerator[AsyncClient, Any]:
    async with AsyncClient(
        app=app, base_url="http://test/api"
    ) as client, LifespanManager(app):
        yield client
