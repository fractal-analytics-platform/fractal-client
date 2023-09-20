from datetime import datetime
from typing import Any
from typing import Optional

from pydantic import BaseModel


__all__ = (
    "_StateBase",
    "StateRead",
)


class _StateBase(BaseModel):
    """
    Base class for `State`

    Attributes:
        id: Primary key
        data: Content of the state
        timestamp: Time stamp of the state
    """

    data: dict[str, Any]
    timestamp: datetime

    def sanitised_dict(self):
        """
        Return `self.dict()` after casting `self.timestamp` to a string
        """
        d = self.dict()
        d["timestamp"] = str(self.timestamp)
        return d


class StateRead(_StateBase):
    id: Optional[int]
