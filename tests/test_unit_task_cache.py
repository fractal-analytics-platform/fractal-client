import pytest
from devtools import debug

from fractal.cmd._aux_task_caching import _get_task_id


TASK_LIST = [
    # Common tasks (with owner=None)
    dict(name="Task A", id=1, version="0.1.1", owner=None),
    dict(name="Task A", id=2, version="1.2.0", owner=None),
    dict(name="Task A", id=3, version="1.3.1", owner=None),
    dict(name="Task B", id=4, version="0.1.1", owner=None),
    dict(name="Task B", id=5, version="1.2.0", owner=None),
    dict(name="Task B", id=6, version="1.3.1", owner=None),
    # UserA tasks
    dict(name="Task A", id=7, version="1.3.1", owner="User1"),
]

NOTSET = "__THIS_VALUE_IS_NOT_SET__"
RAISE = "__THIS_MUST_FAIL__"

CASES = [
    (NOTSET, NOTSET, NOTSET, 1, 1),
    (NOTSET, NOTSET, "User1", NOTSET, 7),
    (NOTSET, NOTSET, NOTSET, 999, RAISE),
    # FIXME: add more cases
]


def test_get_task_id():
    for (name, version, owner, _id, expected_id) in CASES:
        kwargs = {}
        if name != NOTSET:
            kwargs["name"] = name
        if version != NOTSET:
            kwargs["version"] = version
        if _id != NOTSET:
            kwargs["_id"] = _id
        if owner != NOTSET:
            kwargs["owner"] = owner
        debug(kwargs)

        if expected_id == RAISE:
            with pytest.raises(ValueError):
                _get_task_id(TASK_LIST, **kwargs)
        else:
            assert _get_task_id(TASK_LIST, **kwargs) == expected_id
