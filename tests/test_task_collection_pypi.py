import json
from pathlib import Path

import pytest
from devtools import debug

from fractal_server.config import get_settings
from fractal_server.syringe import Inject
from fractal_server.tasks.collection import _init_venv
from fractal_server.tasks.collection import _pip_install
from fractal_server.tasks.collection import _TaskCollectPip
from fractal_server.tasks.collection import create_package_dir_pip
from fractal_server.tasks.collection import download_package
from fractal_server.tasks.collection import inspect_package
from fractal_server.tasks.collection import load_manifest
from fractal_server.tasks.collection import ManifestV1


@pytest.mark.parametrize(
    ("package", "version", "source"),
    [
        ("my_package", None, "pip:my_package"),
        ("my-package", "1.2.3", "pip:my-package==1.2.3"),
    ],
)
def test_unit_source_resolution(package, version, source):
    """
    GIVEN a private task package
    WHEN the source is resolved
    THEN it matches expectations
    """
    tc = _TaskCollectPip(package=package, version=version)
    assert tc.source == source


def test_task_collect_model(dummy_task_package):
    debug(dummy_task_package)
    tc = _TaskCollectPip(package=dummy_task_package.as_posix())

    assert tc.package == "fractal_tasks_dummy"
    assert tc.package_path == dummy_task_package


async def test_init_venv(tmp_path):
    """
    GIVEN a path and a python version
    WHEN _init_venv() is called
    THEN a python venv is initialised at path
    """
    venv_path = tmp_path / "fractal_test"
    venv_path.mkdir(exist_ok=True, parents=True)
    logger_name = "fractal"

    python_bin = await _init_venv(
        path=venv_path, python_version="3.8", logger_name=logger_name
    )

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
    logger_name = "fractal"

    await _init_venv(
        path=venv_path, python_version="3.8", logger_name=logger_name
    )
    location = await _pip_install(
        venv_path=venv_path,
        task_pkg=_TaskCollectPip(package=PACKAGE, version=VERSION),
        logger_name=logger_name,
    )
    debug(location)
    assert PACKAGE in location.as_posix()


async def test_download(tmp_path):
    """
    GIVEN a PyPI package name
    WHEN download_package is called
    THEN the package's wheel is download in the destination directory
    """
    PACKAGE = "fractal-tasks-core"
    task_pkg = _TaskCollectPip(package=PACKAGE)
    pkg = await download_package(task_pkg=task_pkg, dest=tmp_path)
    debug(pkg)
    assert pkg.exists()
    assert "whl" in pkg.as_posix()


async def test_inspect(tmp_path):
    """
    GIVEN the path to a wheel package
    WHEN the inspect package is called on the path of the wheel
    THEN the version number and the manifest are loaded
    """
    PACKAGE = "fractal-tasks-core"
    task_pkg = _TaskCollectPip(package=PACKAGE)
    pkg = await download_package(task_pkg=task_pkg, dest=tmp_path)

    res = inspect_package(pkg)
    debug(res)
    assert "version" in res
    assert isinstance(res["manifest"], ManifestV1)


def test_load_manifest(tmp_path):
    TASK_EXECUTABLE = "task0.py"
    __FRACTAL_MANIFEST__ = dict(
        manifest_version=1,
        task_list=[
            dict(
                name="task0",
                executable=TASK_EXECUTABLE,
                input_type="Any",
                output_type="Any",
                default_args=dict(a=1, b="c"),
                meta=dict(
                    min_memory="2G", requires_gpu=True, version="custom"
                ),
            )
        ],
    )

    package_root = tmp_path / "package"
    package_root.mkdir(exist_ok=True, parents=True)
    manifest_path = package_root / "__FRACTAL_MANIFEST__.json"
    executable_path = package_root / TASK_EXECUTABLE
    python_bin = package_root / "my/custon/python_bin"
    SOURCE = "my:source==123"

    with executable_path.open("w") as f:
        f.write("test_executable")

    with manifest_path.open("w") as f:
        json.dump(__FRACTAL_MANIFEST__, f)

    task_list = load_manifest(
        package_root=package_root,
        python_bin=python_bin,
        source=SOURCE,
    )

    debug(task_list)

    assert len(task_list) == 1
    for t in task_list:
        assert t.source == SOURCE


@pytest.mark.parametrize(
    ("task_pkg", "user", "expected_path"),
    [
        (
            _TaskCollectPip(package="my-package"),
            None,
            Path(".fractal/my-package"),
        ),
        (
            _TaskCollectPip(package="my-package", version="1.2.3"),
            None,
            Path(".fractal/my-package1.2.3"),
        ),
        (
            _TaskCollectPip(package="my-package", version="1.2.3"),
            "genoveffo",
            Path("genoveffo/my-package1.2.3"),
        ),
    ],
)
def test_create_pkg_dir(task_pkg, user, expected_path):
    """
    GIVEN a taks package
    WHEN the directory for installation is created
    THEN the path is the one expected

    NOTE:
        expected_path relative to FRACTAL_ROOT
    """
    settings = Inject(get_settings)
    check = settings.FRACTAL_ROOT / expected_path

    venv_path = create_package_dir_pip(task_pkg=task_pkg, user=user)
    assert venv_path == check
