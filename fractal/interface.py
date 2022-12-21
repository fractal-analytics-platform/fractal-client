from typing import Any
from typing import Dict
from typing import Optional
from typing import Sequence

from rich import print_json
from rich.console import Console


class BaseInterface:
    def __init__(self, retcode: int, data=None):
        self.retcode = retcode
        self.data = data

    def show(self, *args, **kwargs):
        raise NotImplementedError("Implement in subclasses")

    def __repr__(self):
        return f"<Interface (retcode={self.retcode})>"


class PrintInterface(BaseInterface):
    """
    Basic interface that simply prints a string to screen
    """

    def __init__(self, retcode: int, data: Any):
        super().__init__(retcode, data)

    def show(self, *args, **kwargs):
        print(str(self.data))


class RichJsonInterface(BaseInterface):
    """
    Output json using rich.print_json
    """

    def __init__(
        self,
        retcode: int,
        data: Dict[str, Any],
        extra_lines: Optional[Sequence[str]] = None,
    ):
        super().__init__(retcode, data)
        if extra_lines is None:
            self.extra_lines = None
        else:
            self.extra_lines = "\n".join(extra_lines)

    def show(self, *args, **kwargs):
        print_json(data=self.data)
        if self.extra_lines:
            print(self.extra_lines)


class RichConsoleInterface(BaseInterface):
    def __init__(self, retcode: int, data: Any):
        super().__init__(retcode, data)

    def show(self, *args, **kwargs):
        console = Console()
        console.print(self.data)
