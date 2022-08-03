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

from fractal.common.models import SubtaskBase
from fractal.common.models import TaskBase


def flatten(xx):
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

    @property
    def _arguments(self):
        return self.args

    @property
    def callable(self):
        return self.module.partition(":")[2]

    @property
    def import_path(self):
        return self.module.partition(":")[0]


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
        out = self.subtask.default_args.copy()
        out.update(self.args)
        return out

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
                    args=self.default_args,
                )
            ]

    @property
    def _arguments(self):
        if not self._is_atomic:
            raise ValueError("Cannot call _arguments on a non-atomic task")
        return self.default_args

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
        subtask: "Task",
        order: Optional[int] = None,
        args: Optional[Dict[str, Any]] = None,
        commit_and_refresh: bool = True,
    ):
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
