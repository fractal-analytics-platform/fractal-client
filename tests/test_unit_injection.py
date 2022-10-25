import pytest
from devtools import debug

from fractal_server.dependency_injection import Inject


def test_injection():
    debug(Inject)

    INT_VALUE = 5
    Inject.register(int, INT_VALUE)

    debug(Inject._dependencies)

    assert Inject(int) == INT_VALUE

    with pytest.raises(RuntimeError):
        from fractal_server.dependency_injection import _Inject

        _Inject()


def test_injection_override():
    INT_VALUE = 5
    Inject.register(int, 5)
    assert Inject(int) == INT_VALUE

    Inject.register(int, 2 * INT_VALUE)
    assert Inject(int) == INT_VALUE * 2
