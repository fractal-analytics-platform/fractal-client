from typing import Any
from typing import Dict
from typing import Optional

from pydantic import BaseModel


class MockTask(BaseModel):
    name: str
    command: str
    parallelization_level: Optional[str] = None


class MockWorkflowTask(BaseModel):
    order: int = 0
    task: MockTask
    arguments: Dict = {}

    @property
    def is_parallel(self) -> bool:
        return bool(self.task.parallelization_level)

    @property
    def parallelization_level(self) -> Optional[str]:
        return self.task.parallelization_level

    def assemble_args(self, extra: Dict[str, Any] = None):
        """
        Merge of `extra` arguments and `self.arguments`.

        Return
        ------
        full_arsgs (Dict):
            A dictionary consisting of the merge of `extra` and
            self.arguments.
        """
        full_args = {}
        if extra:
            full_args.update(extra)
        full_args.update(self.arguments)
        return full_args
