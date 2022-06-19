import pytest
from fastapi import FastAPI
from typing import AsyncGenerator, Any
from httpx import AsyncClient
from asgi_lifespan import LifespanManager
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
