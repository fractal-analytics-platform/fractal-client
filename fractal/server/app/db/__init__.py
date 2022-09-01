from typing import AsyncGenerator
from typing import Generator

from sqlalchemy import create_engine
from sqlalchemy.ext.asyncio import AsyncSession
from sqlalchemy.ext.asyncio import create_async_engine
from sqlalchemy.orm import Session as DBSyncSession
from sqlalchemy.orm import sessionmaker

from ...config import settings

"""
Losely adapted from https://testdriven.io/blog/fastapi-sqlmodel/#async-sqlmodel
"""


engine = create_async_engine(
    settings.DATABASE_URL, echo=settings.DB_ECHO, future=True
)

engine_sync = create_engine(
    settings.DATABASE_SYNC_URL,
    echo=settings.DB_ECHO,
    future=True,
    connect_args={"check_same_thread": False},
)

async_session_maker = sessionmaker(
    engine, class_=AsyncSession, expire_on_commit=False
)


sync_session_maker = sessionmaker(
    bind=engine_sync, autocommit=False, autoflush=False
)


async def get_db() -> AsyncGenerator[AsyncSession, None]:
    async with async_session_maker() as async_session:
        yield async_session


def get_sync_db() -> Generator[DBSyncSession, None, None]:
    with sync_session_maker() as sync_session:
        yield sync_session
