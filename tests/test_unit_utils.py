from datetime import timezone
from pathlib import Path

import pytest
from devtools import debug

from fractal_server.utils import execute_command
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


async def test_execute_command(tmp_path):
    """
    GIVEN the `pwd` command and a path
    WHEN executed via _execute_command with the path as working directory
    THEN the command returns the full path
    """
    res = await execute_command(cwd=tmp_path, command="pwd")
    assert res.strip() == tmp_path.as_posix()


async def test_execute_command_fail():
    """
    GIVEN a unix command that fails
    WHEN the command is executed via _execute_command
    THEN an error is raised
    """
    with pytest.raises(RuntimeError) as e:
        await execute_command(cwd=Path("/tmp"), command="ls __NOEXIST__")
    debug(e.value)
