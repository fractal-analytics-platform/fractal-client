from typing import AsyncGenerator
from typing import Generator
from warnings import warn

from sqlalchemy import create_engine
from sqlalchemy.ext.asyncio import AsyncSession
from sqlalchemy.ext.asyncio import create_async_engine
from sqlalchemy.orm import Session as DBSyncSession
from sqlalchemy.orm import sessionmaker

from ...config import get_settings
from ...syringe import Inject

"""
Losely adapted from https://testdriven.io/blog/fastapi-sqlmodel/#async-sqlmodel
"""


class DB:
    @classmethod
    def engine_async(cls):
        try:
            return cls._engine_async
        except AttributeError:
            cls.set_db()
            return cls._engine_async

    @classmethod
    def engine_sync(cls):
        try:
            return cls._engine_sync
        except AttributeError:
            cls.set_db()
            return cls._engine_sync

    @classmethod
    def set_db(cls):
        settings = Inject(get_settings)

        if settings.DB_ENGINE == "sqlite":
            warn(
                "SQLite is partially supported but discouraged "
                "in production environment."
                "SQLite offers partial support for ForeignKey "
                "constraints. As such, consistency of the database "
                "cannot be guaranteed."
            )

        cls._engine_async = create_async_engine(
            settings.DATABASE_URL, echo=settings.DB_ECHO, future=True
        )
        cls._engine_sync = create_engine(
            settings.DATABASE_SYNC_URL,
            echo=settings.DB_ECHO,
            future=True,
            connect_args=(
                {"check_same_thread": False}
                if settings.DB_ENGINE == "sqlite"
                else {}
            ),
        )

        cls._async_session_maker = sessionmaker(
            cls._engine_async, class_=AsyncSession, expire_on_commit=False
        )

        cls._sync_session_maker = sessionmaker(
            bind=cls._engine_sync, autocommit=False, autoflush=False
        )

    @classmethod
    async def get_db(cls) -> AsyncGenerator[AsyncSession, None]:
        try:
            session_maker = cls._async_session_maker()
        except AttributeError:
            cls.set_db()
            session_maker = cls._async_session_maker()
        async with session_maker as async_session:
            yield async_session

    @classmethod
    def get_sync_db(cls) -> Generator[DBSyncSession, None, None]:
        with cls._sync_session_maker() as sync_session:
            yield sync_session


get_db = DB.get_db
get_sync_db = DB.get_sync_db
