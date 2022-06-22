from sqlalchemy.ext.asyncio import AsyncSession
from sqlalchemy.ext.asyncio import create_async_engine
from sqlalchemy.orm import sessionmaker

"""
Losely adapted from https://testdriven.io/blog/fastapi-sqlmodel/#async-sqlmodel
"""


DATABASE_URL = "sqlite+aiosqlite://"


engine = create_async_engine(DATABASE_URL, echo=True, future=True)
async_session_maker = sessionmaker(
    engine, class_=AsyncSession, expire_on_commit=False
)


async def init_db():
    pass


async def get_db() -> AsyncSession:
    async with async_session_maker() as session:
        yield session
