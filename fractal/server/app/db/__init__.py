from sqlalchemy.ext.asyncio import AsyncSession
from sqlalchemy.ext.asyncio import create_async_engine
from sqlalchemy.orm import sessionmaker

from ...config import settings

"""
Losely adapted from https://testdriven.io/blog/fastapi-sqlmodel/#async-sqlmodel
"""


engine = create_async_engine(
    settings.DATABASE_URL, echo=settings.DB_ECHO, future=True
)

async_session_maker = sessionmaker(
    engine, class_=AsyncSession, expire_on_commit=False
)


async def get_db() -> AsyncSession:
    async with async_session_maker() as session:
        yield session
