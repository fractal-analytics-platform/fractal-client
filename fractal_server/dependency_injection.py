from typing import Any
from typing import Dict


_instance_count = 0


class _Inject:
    _dependencies: Dict[Any, Any] = {}

    def __init__(self):
        global _instance_count
        if _instance_count == 1:
            raise RuntimeError("You must only instance this class once")
        _instance_count += 1

    @classmethod
    def __call__(cls, Type):
        return cls._dependencies[Type]

    @classmethod
    def register(cls, Type, value) -> None:
        cls._dependencies[Type] = value


Inject = _Inject()
