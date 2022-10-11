from collections.abc import MutableMapping
from typing import Any


def popget(d: MutableMapping, key: str, default: Any = None) -> Any:
    """
    Pop and return mapping item if possible or return default
    """
    try:
        return d.pop(key)
    except KeyError:
        return default
