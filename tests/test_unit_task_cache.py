import pytest

from fractal.cmd._aux_task_caching import _get_matching_tasks
from fractal.cmd._aux_task_caching import _search_in_task_list
from fractal.cmd._aux_task_caching import FractalCacheError


def test_get_matching_tasks(clear_task_cache):
    """Test all possible cases for function `_get_matching_tasks`"""

    a = dict(name="bob", id=1, version="0.1.1")
    b = dict(name="bob", id=2, version="1.2.0")
    c = dict(name="bob", id=3, version="1.3.1")
    d = dict(name="bar", id=4, version="0.0.0")
    e = dict(name="bar", id=5, version=None)

    TASK_LIST = [a, b, c, d, e]

    res = _get_matching_tasks(TASK_LIST, name="alice")
    assert res == []

    res = _get_matching_tasks(TASK_LIST, name="bob")
    assert res == [a, b, c]

    res = _get_matching_tasks(TASK_LIST, name="bob", version="1.2.0")
    assert res == [b]

    res = _get_matching_tasks(TASK_LIST, name="bar")
    assert res == [d, e]

    res = _get_matching_tasks(TASK_LIST, name="bar", version="3.1.4")
    assert res == []


async def test_search_in_task_list(clear_task_cache):
    """Test all possible cases for function `_search_in_task_list`"""

    TASK_LIST = [
        dict(name="dummy1", id=101, version="1.0.1", source="a"),
        dict(name="dummy2", id=201, version=None, source="b"),
        dict(name="dummy2", id=202, version="2.0.0", source="c"),
        dict(name="dummy3", id=301, version="3.0.0", source="d"),
        dict(name="dummy3", id=302, version="3.1.4", source="e"),
        dict(name="dummy4", id=401, version="4.0.0", source="f"),
        dict(name="dummy4", id=402, version="4.1.1", source="g"),
        dict(name="dummy4", id=401, version="4.1.1", source="h"),
    ]

    # TEST zero matching

    # case 1
    with pytest.raises(FractalCacheError) as err:
        res = _search_in_task_list(task_list=TASK_LIST, name="dummy0")
    print(err.value.args[0])
    assert 'There is no task with name "dummy0"' in err.value.args[0]  # noqa

    # case 2
    with pytest.raises(FractalCacheError) as err:
        res = _search_in_task_list(
            task_list=TASK_LIST, name="dummy1", version="3.1.4"
        )
    print(err.value.args[0])
    assert (
        'There is no task with (name, version)=("dummy1", 3.1.4)'
        in err.value.args[0]
    )  # noqa

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
    print(err.value.args[0])
    assert "Cannot determine the latest version" in err.value.args[0]
    # case 2
    res = _search_in_task_list(task_list=TASK_LIST, name="dummy3")
    assert res == 302
    # case 3
    with pytest.raises(FractalCacheError) as err:
        res = _search_in_task_list(task_list=TASK_LIST, name="dummy4")
    print(err.value.args[0])
    assert "Multiple tasks with latest version (4.1.1)" in err.value.args[0]
    print(err.value.args[0])
    # case 4
    with pytest.raises(FractalCacheError) as err:
        res = _search_in_task_list(
            task_list=TASK_LIST, name="dummy4", version="4.1.1"
        )
    print(err.value.args[0])
    assert "Multiple tasks with version 4.1.1" in err.value.args[0]
