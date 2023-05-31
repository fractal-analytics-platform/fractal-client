import pytest
from devtools import debug

from fractal.cmd._aux_task_caching import _get_matching_tasks
from fractal.cmd._aux_task_caching import _search_in_task_list
from fractal.cmd._aux_task_caching import FractalCacheError


def test_get_task_id(clear_task_cache):
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
        (NOTSET, NOTSET, 1, 1),
        (NOTSET, NOTSET, NOTSET, 7),
        (NOTSET, NOTSET, 999, RAISE),
        # FIXME: add more cases
    ]
    for (name, version, _id, expected_id) in CASES:
        kwargs = {}
        if name != NOTSET:
            kwargs["name"] = name
        if version != NOTSET:
            kwargs["version"] = version
        if _id != NOTSET:
            kwargs["_id"] = _id
        debug(kwargs)

        if expected_id == RAISE:
            with pytest.raises(ValueError):
                _get_matching_tasks(TASK_LIST, **kwargs)
        else:
            assert _get_matching_tasks(TASK_LIST, **kwargs) == expected_id


async def test_search_in_task_list(clear_task_cache):

    TASK_LIST = [
        dict(name="dummy1", id=101, version="1.0.1"),
        dict(name="dummy2", id=201, version=None),
        dict(name="dummy2", id=202, version="2.0.0"),
        dict(name="dummy3", id=301, version="3.0.0"),
        dict(name="dummy3", id=302, version="3.1.4"),
        dict(name="dummy4", id=401, version="4.0.0"),
        dict(name="dummy4", id=402, version="4.1.1"),
        dict(name="dummy4", id=401, version="4.1.1"),
    ]

    # TEST zero matching
    # case 1
    with pytest.raises(FractalCacheError) as err:
        res = _search_in_task_list(task_list=TASK_LIST, name="dummy0")
    assert err.value.args[0] == (
        "There is no task with (name, version)=(dummy0,None) in the cache"
    )
    # case 2
    with pytest.raises(FractalCacheError) as err:
        res = _search_in_task_list(
            task_list=TASK_LIST, name="dummy1", version="3.1.4"
        )
    assert err.value.args[0] == (
        "There is no task with (name, version)=(dummy1,3.1.4) in the cache"
    )
    # TEST one matching
    # case 1
    res = _search_in_task_list(task_list=TASK_LIST, name="dummy1")
    assert res == 101
    # case 2
    res = _search_in_task_list(
        task_list=TASK_LIST, name="dummy1", version="1.0.1"
    )
    assert res == 101

    # TEST multiple matching
    # case 1
    with pytest.raises(FractalCacheError) as err:
        res = _search_in_task_list(task_list=TASK_LIST, name="dummy2")
    assert (
        "Cannot determine the maximum version in this list "
        "(there are one or more None)\n"
    ) in err.value.args[0]
    # case 2
    res = _search_in_task_list(task_list=TASK_LIST, name="dummy3")
    assert res == 302
    # case 3
    with pytest.raises(FractalCacheError) as err:
        res = _search_in_task_list(task_list=TASK_LIST, name="dummy4")
    assert "Multiple tasks with maximum version 4.1.1\n" in err.value.args[0]
    # case 4
    with pytest.raises(FractalCacheError) as err:
        res = _search_in_task_list(
            task_list=TASK_LIST, name="dummy4", version="4.1.1"
        )
    assert "Multiple tasks with version 4.1.1\n" in err.value.args[0]
