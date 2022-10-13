from collections.abc import MutableMapping
from datetime import datetime
from datetime import timezone
from typing import Any


def get_timestamp() -> datetime:
    return datetime.now(tz=timezone.utc)


def popget(d: MutableMapping, key: str, default: Any = None) -> Any:
    """
    Pop and return mapping item if possible or return default
    """
    try:
        return d.pop(key)
    except KeyError:
        return default
