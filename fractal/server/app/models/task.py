from typing import Any
from typing import Dict
from typing import List
from typing import Optional

from pydantic import BaseModel
from sqlalchemy import Column
from sqlalchemy.ext.asyncio import AsyncSession
from sqlalchemy.ext.orderinglist import ordering_list
from sqlalchemy.types import JSON
from sqlmodel import Field
from sqlmodel import Relationship

from ...utils import popget
from fractal.common.models import SubtaskBase
from fractal.common.models import TaskBase


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


class Subtask(SubtaskBase, table=True):  # type: ignore
    class Config:
        arbitrary_types_allowed = True
        fields = {"parent": {"exclude": True}}

    parent_task_id: Optional[int] = Field(
        default=None, foreign_key="task.id", primary_key=True
    )
    subtask_id: Optional[int] = Field(
        default=None, foreign_key="task.id", primary_key=True
    )
    parent: "Task" = Relationship(
        sa_relationship_kwargs=dict(
            primaryjoin="Subtask.parent_task_id==Task.id"
        )
    )
    subtask: "Task" = Relationship(
        sa_relationship_kwargs=dict(primaryjoin="Subtask.subtask_id==Task.id")
    )
    order: Optional[int]
    args: Dict[str, Any] = Field(sa_column=Column(JSON), default={})

    @property
    def _is_atomic(self):
        """
        A subtask is atomic iff the child task is atomic
        """
        return self.subtask._is_atomic

    @property
    def _arguments(self):
        """
        Override default arguments and strip specific arguments (executor and
        parallelization_level)
        """
        out = self.subtask.default_args.copy()
        out.update(self.args)
        popget(out, "executor")
        popget(out, "parallelization_level")
        return out

    @property
    def executor(self) -> Optional[str]:
        return self.args.get("executor") or self.subtask.executor

    @property
    def parallelization_level(self) -> Optional[str]:
        return (
            self.args.get("parallelization_level")
            or self.subtask.parallelization_level
        )

    @property
    def import_path(self):
        return self.subtask.import_path

    @property
    def callable(self):
        return self.subtask.callable

    def preprocess(self):
        if self._is_atomic:
            return PreprocessedTask(
                name=self.subtask.name,
                module=self.subtask.module,
                args=self._arguments,
                executor=self.executor,
                parallelization_level=self.parallelization_level,
            )
        else:
            return [st.preprocess() for st in self.subtask.subtask_list]

    @property
    def name(self):
        return self.subtask.name


class Task(TaskBase, table=True):  # type: ignore
    id: Optional[int] = Field(default=None, primary_key=True)
    default_args: Dict[str, Any] = Field(sa_column=Column(JSON), default={})
    subtask_list: List["Subtask"] = Relationship(
        sa_relationship_kwargs=dict(
            primaryjoin="Task.id==Subtask.parent_task_id",
            lazy="selectin",
            order_by="Subtask.order",
            collection_class=ordering_list("order"),
        ),
    )

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

    async def add_subtask(
        self,
        db: AsyncSession,
        subtask_id: Optional[int] = None,
        subtask: Optional["Task"] = None,
        order: Optional[int] = None,
        args: Optional[Dict[str, Any]] = None,
        commit_and_refresh: bool = True,
    ):
        if subtask is None:
            subtask = await db.get(Task, subtask_id)

        if not args:
            args = dict()
        st = Subtask(parent=self, subtask=subtask, args=args)

        if order is None:
            self.subtask_list.append(st)
        else:
            self.subtask_list.insert(order, st)
        db.add_all([self, st])
        if commit_and_refresh:
            await db.commit()
            await db.refresh(self)
