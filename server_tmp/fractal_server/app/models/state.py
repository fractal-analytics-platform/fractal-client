from datetime import datetime
from typing import Any
from typing import Dict
from typing import Optional

from sqlalchemy import Column
from sqlalchemy.types import DateTime
from sqlalchemy.types import JSON
from sqlmodel import Field
from sqlmodel import SQLModel

from ...utils import get_timestamp


class State(SQLModel, table=True):
    """
    Store arbitrary data in the database

    This table is just a state interchange that allows the system to store
    arbitrary data for later retrieval. This is particuarly important for long
    background tasks, in which it is not possible to return a meaningful
    response to the client within a single request lifespan.
    """

    id: Optional[int] = Field(default=None, primary_key=True)
    data: Dict[str, Any] = Field(sa_column=Column(JSON), default={})
    timestamp: datetime = Field(
        default_factory=get_timestamp,
        nullable=False,
        sa_column=Column(DateTime(timezone=True)),
    )

    class Config:
        arbitrary_types_allowed = True
