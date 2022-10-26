from typing import Any
from typing import Dict
from typing import Optional

from sqlalchemy import Column
from sqlalchemy.types import JSON
from sqlmodel import Field

from ..schemas.task import _TaskBase


class Task(_TaskBase, table=True):  # type: ignore
    id: Optional[int] = Field(default=None, primary_key=True)
    project_id: Optional[int] = Field(foreign_key="project.id")
    default_args: Dict[str, Any] = Field(sa_column=Column(JSON), default={})

    @property
    def parallelization_level(self) -> Optional[str]:
        try:
            return self.default_args["parallelization_level"]
        except KeyError:
            return None

    @property
    def is_parallel(self):
        return bool(self.parallelization_level)

    @property
    def executor(self) -> Optional[str]:
        try:
            return self.default_args["executor"]
        except KeyError:
            return None
