import json
from typing import Any


class Interface:
    def __init__(self, retcode: int, data: Any) -> None:
        self.retcode = retcode
        self.data = data

    def show(self, *args, **kwargs):
        if isinstance(self.data, (dict, list)):
            print(json.dumps(self.data, indent=2, sort_keys=True))
        else:
            print(self.data)
