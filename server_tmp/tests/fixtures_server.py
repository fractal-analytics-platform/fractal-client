"""
Copyright 2022 (C) Friedrich Miescher Institute for Biomedical Research and
University of Zurich

Original authors:
Jacopo Nespolo <jacopo.nespolo@exact-lab.it>

This file is part of Fractal and was originally developed by eXact lab S.r.l.
<exact-lab.it> under contract with Liberali Lab from the Friedrich Miescher
Institute for Biomedical Research and Pelkmans Lab from the University of
Zurich.
"""
import asyncio
import logging
from dataclasses import dataclass
from dataclasses import field
from pathlib import Path
from typing import Any
from typing import AsyncGenerator
from typing import List
from typing import Optional
from uuid import uuid4

import pytest
from asgi_lifespan import LifespanManager
from fastapi import FastAPI
from httpx import AsyncClient
from sqlalchemy.ext.asyncio import AsyncSession

from fractal_server.config import get_settings
from fractal_server.config import Settings
from fractal_server.syringe import Inject


DB_ENGINE = "postgres"


def get_patched_settings(temp_path: Path):
    settings = Settings()
    settings.JWT_SECRET_KEY = "secret_key"
    settings.DEPLOYMENT_TYPE = "development"

    settings.DB_ENGINE = DB_ENGINE
    if DB_ENGINE == "sqlite":
        settings.SQLITE_PATH = (
            f"{temp_path.as_posix()}/_test.db?mode=memory&cache=shared"
        )
    elif DB_ENGINE == "postgres":
        settings.DB_ENGINE = "postgres"
        settings.POSTGRES_USER = "postgres"
        settings.POSTGRES_PASSWORD = "postgres"
        settings.POSTGRES_DB = "fractal"
    else:
        raise ValueError

    settings.FRACTAL_ROOT = temp_path
    settings.RUNNER_ROOT_DIR = temp_path / "artifacts"
    settings.FRACTAL_LOGGING_LEVEL = logging.DEBUG
    return settings


@pytest.fixture(scope="session", autouse=True)
def override_settings(tmp_path_factory):
    tmp_path = tmp_path_factory.mktemp("fractal_root")

    def _get_settings():
        return get_patched_settings(tmp_path)

    Inject.override(get_settings, _get_settings)
    try:
        yield
    finally:
        Inject.pop(get_settings)


@pytest.fixture
def unset_deployment_type():
    """
    Temporarily override the seetings with a version that would fail
    `settings.check()`

    Afterwards, restore any previous injection, if any.
    """
    settings = Settings()
    settings.DEPLOYMENT_TYPE = None

    def _get_settings():
        return settings

    try:
        previous = Inject.pop(get_settings)
    except RuntimeError:
        previous = None

    Inject.override(get_settings, _get_settings)
    try:
        yield
    finally:
        Inject.pop(get_settings)
        if previous:
            Inject.override(get_settings, previous)


@pytest.fixture(scope="session")
def event_loop():
    _event_loop = asyncio.new_event_loop()
    yield _event_loop


@pytest.fixture
async def db_create_tables(override_settings):
    from fractal_server.app.db import DB
    from fractal_server.app.models import SQLModel

    engine = DB.engine_sync()
    metadata = SQLModel.metadata
    metadata.create_all(engine)
    yield

    metadata.drop_all(engine)


@pytest.fixture
async def db(db_create_tables):
    from fractal_server.app.db import get_db

    async for session in get_db():
        yield session


@pytest.fixture
async def db_sync(db_create_tables):
    from fractal_server.app.db import get_sync_db

    for session in get_sync_db():
        yield session


@pytest.fixture
async def app(override_settings) -> AsyncGenerator[FastAPI, Any]:
    app = FastAPI()
    yield app


@pytest.fixture
async def register_routers(app, override_settings):
    from fractal_server.main import collect_routers

    collect_routers(app)


@pytest.fixture
async def client(
    app: FastAPI, register_routers, db
) -> AsyncGenerator[AsyncClient, Any]:
    async with AsyncClient(
        app=app, base_url="http://test"
    ) as client, LifespanManager(app):
        yield client


@pytest.fixture
async def MockCurrentUser(app, db):
    from fractal_server.app.security import current_active_user
    from fractal_server.app.security import User

    @dataclass
    class _MockCurrentUser:
        """
        Context managed user override
        """

        name: str = "User Name"
        scopes: Optional[List[str]] = field(
            default_factory=lambda: ["project"]
        )
        email: Optional[str] = field(
            default_factory=lambda: f"{uuid4()}@exact-lab.it"
        )
        persist: Optional[bool] = False

        def _create_user(self):
            self.user = User(
                name=self.name,
                email=self.email,
                hashed_password="fake_hashed_password",
                slurm_user="slurm_user",
            )

        def current_active_user_override(self):
            def __current_active_user_override():
                return self.user

            return __current_active_user_override

        async def __aenter__(self):
            self._create_user()

            if self.persist:
                db.add(self.user)
                await db.commit()
                await db.refresh(self.user)
                # Removing object from test db session, so that we can operate
                # on user from other sessions
                db.expunge(self.user)
            self.previous_user = app.dependency_overrides.get(
                current_active_user, None
            )
            app.dependency_overrides[
                current_active_user
            ] = self.current_active_user_override()
            return self.user

        async def __aexit__(self, *args, **kwargs):
            if self.previous_user:
                app.dependency_overrides[
                    current_active_user
                ] = self.previous_user()

    return _MockCurrentUser


@pytest.fixture
async def project_factory(db):
    """
    Factory that adds a project to the database
    """
    from fractal_server.app.models import Project

    async def __project_factory(user, **kwargs):
        defaults = dict(
            name="project",
            project_dir="/tmp/",
            slug="slug",
        )
        defaults.update(kwargs)
        project = Project(**defaults)
        project.user_member_list.append(user)
        db.add(project)
        await db.commit()
        await db.refresh(project)
        return project

    return __project_factory


@pytest.fixture
async def dataset_factory(db):
    from fractal_server.app.models import Project, Dataset

    async def __dataset_factory(project: Project, **kwargs):
        defaults = dict(name="test dataset")
        defaults.update(kwargs)
        project.dataset_list.append(Dataset(**defaults))
        db.add(project)
        await db.commit()
        await db.refresh(project)
        return project.dataset_list[-1]

    return __dataset_factory


@pytest.fixture
async def resource_factory(db, testdata_path):
    from fractal_server.app.models import Dataset, Resource

    async def __resource_factory(dataset: Dataset, **kwargs):
        """
        Add a new resorce to dataset
        """
        defaults = dict(
            path=(testdata_path / "png").as_posix(), glob_pattern="*.png"
        )
        defaults.update(kwargs)
        resource = Resource(dataset_id=dataset.id, **defaults)
        db.add(resource)
        await db.commit()
        await db.refresh(dataset)
        return dataset.resource_list[-1]

    return __resource_factory


@pytest.fixture
async def task_factory(db: AsyncSession):
    """
    Insert task in db
    """
    from fractal_server.app.models import Task

    async def __task_factory(db: AsyncSession = db, index: int = 0, **kwargs):
        defaults = dict(
            name=f"task{index}",
            input_type="zarr",
            output_type="zarr",
            command="cmd",
            source="source",
        )
        args = dict(**defaults)
        args.update(kwargs)
        t = Task(**args)
        db.add(t)
        await db.commit()
        await db.refresh(t)
        return t

    return __task_factory
