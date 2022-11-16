from datetime import datetime
from typing import Any
from typing import Dict
from typing import Optional

from sqlmodel import SQLModel

__all__ = (
    "_StateBase",
    "StateRead",
)


class _StateBase(SQLModel):
    id: Optional[int]
    data: Dict[str, Any]
    timestamp: datetime

    class Config:
        arbitrary_types_allowed = True

    def sanitised_dict(self):
        d = self.dict()
        d["timestamp"] = str(self.timestamp)
        return d


class StateRead(_StateBase):
    pass
