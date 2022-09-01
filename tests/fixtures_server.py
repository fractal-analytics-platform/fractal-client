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
from dataclasses import dataclass
from dataclasses import field
from typing import Any
from typing import AsyncGenerator
from typing import List
from typing import Optional
from uuid import uuid4

import pytest
from asgi_lifespan import LifespanManager
from fastapi import FastAPI
from httpx import AsyncClient
from sqlalchemy.ext.asyncio import AsyncEngine
from sqlalchemy.ext.asyncio import AsyncSession
from sqlalchemy.orm import sessionmaker


def override_enironment(testdata_path):
    from os import environ

    environ["JWT_SECRET_KEY"] = "secret_key"
    environ["DEPLOYMENT_TYPE"] = "development"
    environ["DATA_DIR_ROOT"] = testdata_path.as_posix()

    environ["DB_ENGINE"] = "sqlite"
    # Shared in memory database,
    # c.f., https://stackoverflow.com/a/38089822/283972
    environ["SQLITE_PATH"] = "file:cachedb?mode=memory&cache=shared"

    from fractal.server.config import settings

    return settings


@pytest.fixture(scope="session")
def event_loop():
    _event_loop = asyncio.new_event_loop()
    yield _event_loop


@pytest.fixture(autouse=True, scope="session")
def patch_settings(testdata_path):
    return override_enironment(testdata_path)


@pytest.fixture(scope="session")
async def db_engine(patch_settings) -> AsyncGenerator[AsyncEngine, None]:
    from fractal.server.app.db import engine

    yield engine


@pytest.fixture(scope="session")
def db_sync_engine(patch_settings):
    from fractal.server.app.db import engine_sync

    yield engine_sync


@pytest.fixture
async def db_session_maker(
    db_engine, app
) -> AsyncGenerator[AsyncSession, None]:
    import fractal.server.app.models  # noqa F401 make sure models are imported
    from sqlmodel import SQLModel

    async with db_engine.begin() as conn:
        await conn.run_sync(SQLModel.metadata.create_all)
        async_session_maker = sessionmaker(
            db_engine, class_=AsyncSession, expire_on_commit=False
        )

        async def _get_db():
            async with async_session_maker() as session:
                yield session

        from fractal.server.app.db import get_db

        app.dependency_overrides[get_db] = _get_db

        yield async_session_maker

        await conn.run_sync(SQLModel.metadata.drop_all)


@pytest.fixture
def db_sync_session_maker(db_sync_engine, app):
    from fractal.server.app.db import get_sync_db
    from sqlmodel import Session

    def _get_sync_db():
        with Session(db_sync_engine) as session:
            yield session

    app.dependency_overrides[get_sync_db] = _get_sync_db


@pytest.fixture
async def db(db_session_maker):
    async with db_session_maker() as session:
        yield session


@pytest.fixture()
def db_sync(db_sync_engine):
    from sqlmodel import Session

    with Session(db_sync_engine) as session:
        from devtools import debug

        debug(f"yielding session {id(session)}")
        yield session


@pytest.fixture
async def app(patch_settings) -> AsyncGenerator[FastAPI, Any]:
    app = FastAPI()
    yield app


@pytest.fixture
async def register_routers(app):
    from fractal.server import collect_routers

    collect_routers(app)


@pytest.fixture
async def collect_tasks(db):
    from fractal.server.app.api.v1.task import collect_tasks_headless

    await collect_tasks_headless()


@pytest.fixture
async def client(
    app: FastAPI, register_routers, db, db_sync
) -> AsyncGenerator[AsyncClient, Any]:
    async with AsyncClient(
        app=app, base_url="http://test"
    ) as client, LifespanManager(app):
        yield client


@pytest.fixture
async def MockCurrentUser(app, db):
    from fractal.server.app.security import current_active_user
    from fractal.server.app.security import User

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
    from fractal.server.app.models import Project

    async def __project_factory(user, **kwargs):
        defaults = dict(
            name="project",
            project_dir="/tmp/",
            user_owner_id=user.id,
            slug="slug",
        )
        defaults.update(kwargs)
        project = Project(**defaults)
        db.add(project)
        await db.commit()
        await db.refresh(project)
        return project

    return __project_factory


@pytest.fixture
async def dataset_factory(db):
    from fractal.server.app.models import Project, Dataset

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
    from fractal.server.app.models import Dataset, Resource

    async def __resource_factory(dataset: Dataset, **kwargs):
        defaults = dict(
            path=(testdata_path / "png").as_posix(), glob_pattern="*.png"
        )
        defaults.update(kwargs)
        dataset.resource_list.append(Resource(**defaults))
        db.add(dataset)
        await db.commit()
        await db.refresh(dataset)
        return dataset.resource_list[-1]

    return __resource_factory


@pytest.fixture
async def task_factory(db: AsyncSession):
    """
    Insert task in db
    """
    from fractal.server.app.models import Task

    async def __task_factory(db: AsyncSession = db, index: int = 0, **kwargs):
        defaults = dict(
            name=f"task{index}",
            resource_type="task",
            module=f"task{index}",
            input_type="zarr",
            output_type="zarr",
            subtask_list=[],
        )
        args = dict(**defaults)
        args.update(kwargs)
        subtask_list = args.pop("subtask_list")
        t = Task(**args)
        db.add(t)
        for st in subtask_list:
            await t.add_subtask(db=db, subtask=st, commit_and_refresh=False)
        await db.commit()
        await db.refresh(t)
        return t

    return __task_factory
