from datetime import timezone

from devtools import debug

from fractal_server.utils import get_timestamp


def test_timestamp():
    """
    GIVEN a function that provides a timestamp
    WHEN called
    THEN the timestamp is timezone aware and the timezone is set to UTC
    """
    ts = get_timestamp()
    debug(ts)
    assert ts
    assert ts.tzinfo is timezone.utc
