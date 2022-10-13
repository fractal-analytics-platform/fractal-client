from typing import Any
from typing import Dict
from typing import Optional

from pydantic import BaseModel
from sqlalchemy import Column
from sqlalchemy.types import JSON
from sqlmodel import Field

from ..schemas import TaskBase
from .models_utils import popget


def flatten(xx):
    """
    Flatten an arbitrarily nested sequence of sequences
    """
    for x in xx:
        if isinstance(x, PreprocessedTask):
            yield x
        else:
            yield from x


class PreprocessedTask(BaseModel):
    name: str
    module: str
    args: Dict[str, Any]
    save_intermediate_result: bool = False
    executor: Optional[str] = None
    parallelization_level: Optional[str] = None

    @property
    def callable(self):
        return self.module.partition(":")[2]

    @property
    def import_path(self):
        return self.module.partition(":")[0]

    @property
    def _arguments(self):
        return self.args


class Task(TaskBase, table=True):  # type: ignore
    id: Optional[int] = Field(default=None, primary_key=True)
    project_id: Optional[int] = Field(foreign_key="project.id")
    default_args: Dict[str, Any] = Field(sa_column=Column(JSON), default={})

    def preprocess(self):
        if not self._is_atomic:
            return flatten(st.preprocess() for st in self.subtask_list)
        else:
            return [
                PreprocessedTask(
                    name=self.name,
                    module=self.module,
                    args=self._arguments,
                    executor=self.executor,
                    parallelization_level=self.parallelization_level,
                )
            ]

    @property
    def _arguments(self):
        if not self._is_atomic:
            raise ValueError("Cannot call _arguments on a non-atomic task")
        out = self.default_args.copy()
        popget(out, "executor")
        popget(out, "parallelization_level")
        return out

    @property
    def executor(self) -> Optional[str]:
        try:
            return self.default_args["executor"]
        except KeyError:
            return None

    @property
    def parallelization_level(self) -> Optional[str]:
        try:
            return self.default_args["parallelization_level"]
        except KeyError:
            return None

    @property
    def module(self):
        return self.subtask.module

    @property
    def _is_atomic(self) -> bool:
        """
        A task is atomic iff it does not contain subtasks
        """
        return not self.subtask_list

    @property
    def callable(self):
        return self.module.partition(":")[2]

    @property
    def import_path(self):
        return self.module.partition(":")[0]
