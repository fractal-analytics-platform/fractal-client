from typing import Any
from typing import Dict

from rich import print_json
from rich.console import Console


class BaseInterface:
    def __init__(self, retcode: int):
        self.retcode = retcode

    def show(self, *args, **kwargs):
        raise NotImplementedError("Implement in subclasses")

    def __repr__(self):
        return f"<Interface (retcode={self.retcode})>"


class PrintInterface(BaseInterface):
    """
    Basic interface that simply prints a string to screen
    """

    def __init__(self, retcode: int, output: str):
        super().__init__(retcode)
        self.output = output

    def show(self, *args, **kwargs):
        print(self.output)


class RichJsonInterface(BaseInterface):
    """
    Output json using rich.print_json
    """

    def __init__(self, retcode: int, data: Dict[str, Any]):
        super().__init__(retcode)
        self.data = data

    def show(self, *args, **kwargs):
        print_json(data=self.data)


class RichConsoleInterface(BaseInterface):
    def __init__(self, retcode: int, objects: Any):
        super().__init__(retcode)
        self.objects = objects

    def show(self, *args, **kwargs):
        console = Console()
        console.print(self.objects)
