import pytest
from devtools import debug

from fractal_server.syringe import Inject


THE_ANSWER = 42


def foo():
    debug("Calling foo()")
    return THE_ANSWER


def test_injection():
    assert Inject(foo) == THE_ANSWER


def test_singleton():
    with pytest.raises(RuntimeError):
        from fractal_server.syringe import _Inject

        _Inject()


def test_injection_override():
    def bar(answer=Inject(foo)):
        debug("calling bar()")
        return answer

    assert bar() == THE_ANSWER

    def oof():
        debug("calling oof()")
        return 24

    Inject.override(foo, oof)

    def baz(answer=Inject(foo)):
        debug("calling baz()")
        return answer

    assert baz() == oof()
