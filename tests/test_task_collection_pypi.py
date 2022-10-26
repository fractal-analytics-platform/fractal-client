from pathlib import Path

import pytest
from devtools import debug

from fractal_server.tasks.collection import _execute_command
from fractal_server.tasks.collection import _init_venv
from fractal_server.tasks.collection import _pip_install


async def test_execute_command(tmp_path):
    """
    GIVEN the `pwd` command and a path
    WHEN executed via _execute_command with the path as working directory
    THEN the command returns the full path
    """
    res = await _execute_command(cwd=tmp_path, command="pwd")
    assert res.strip() == tmp_path.as_posix()


async def test_execute_command_fail():
    """
    GIVEN a unix command that fails
    WHEN the command is executed via _execute_command
    THEN an error is raised
    """
    with pytest.raises(RuntimeError) as e:
        await _execute_command(cwd=Path("/tmp"), command="ls __NOEXIST__")
    debug(e.value)


async def test_init_venv(tmp_path):
    """
    GIVEN a path and a python version
    WHEN _init_venv() is called
    THEN a python venv is initialised at path
    """
    venv_path = tmp_path / "fractal_test"
    venv_path.mkdir(exist_ok=True, parents=True)

    python_bin = await _init_venv(path=venv_path, python_version="3.8")

    assert venv_path.exists()
    assert (venv_path / "venv").exists()
    assert (venv_path / "venv/bin/python").exists()
    assert (venv_path / "venv/bin/pip").exists()
    assert python_bin.exists()
    assert python_bin == venv_path / "venv/bin/python"


async def test_pip_install(tmp_path):
    """
    GIVEN a package name and version and a path with a venv
    WHEN _pip_install() is called
    THEN the package is installed in the venv and the package installation
         location is returned
    """
    PACKAGE = "devtools"
    VERSION = "0.8.0"
    venv_path = tmp_path / "fractal_test" / f"{PACKAGE}{VERSION}"
    venv_path.mkdir(exist_ok=True, parents=True)

    await _init_venv(path=venv_path, python_version="3.8")
    location = await _pip_install(
        venv_path=venv_path, package=PACKAGE, version=VERSION
    )
    debug(location)
    assert PACKAGE in location.as_posix()
