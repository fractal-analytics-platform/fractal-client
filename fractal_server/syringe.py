"""
Copyright 2022 (C) eXact lab S.r.l.

Original author:
Jacopo Nespolo <jacopo.nespolo@exact-lab.it>

This file is part of `Syringe`, a simple Python module for dependency
injection. Redistribution and use must comply with the 3-clause BSD License
(see `LICENSE` file).
"""
from typing import Any
from typing import Callable
from typing import Dict
from typing import TypeVar


T = TypeVar("T")
_instance_count = 0


class _Inject:
    _dependencies: Dict[Any, Any] = {}

    def __init__(self):
        global _instance_count
        if _instance_count == 1:
            raise RuntimeError("You must only instance this class once")
        _instance_count += 1

    @classmethod
    def __call__(cls, _callable: Callable[..., T]) -> T:
        try:
            return cls._dependencies[_callable]()
        except KeyError:
            return _callable()

    @classmethod
    def pop(cls, _callable: Callable[..., T]) -> T:
        try:
            return cls._dependencies.pop(_callable)
        except KeyError:
            raise RuntimeError(f"No dependency override for {_callable}")

    @classmethod
    def override(
        cls, _callable: Callable[..., T], value: Callable[..., T]
    ) -> None:
        cls._dependencies[_callable] = value


# NOTE: This is a singleton instance
Inject = _Inject()
