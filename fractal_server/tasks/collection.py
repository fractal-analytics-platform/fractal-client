# Copyright 2022 (C) Friedrich Miescher Institute for Biomedical Research and
# University of Zurich
#
# Original authors:
# Jacopo Nespolo <jacopo.nespolo@exact-lab.it>
#
# This file is part of Fractal and was originally developed by eXact lab S.r.l.
# <exact-lab.it> under contract with Liberali Lab from the Friedrich Miescher
# Institute for Biomedical Research and Pelkmans Lab from the University of
# Zurich.
import json
import logging
from io import IOBase
from pathlib import Path
from typing import List
from typing import Optional
from typing import Tuple
from typing import Union
from zipfile import ZipFile

from pydantic import root_validator

from ..app.schemas import ManifestV1
from ..app.schemas import TaskCollectPip
from ..app.schemas import TaskCreate
from ..config import get_settings
from ..syringe import Inject
from ..utils import close_logger
from ..utils import execute_command
from ..utils import set_logger


class TaskCollectionError(RuntimeError):
    pass


def get_collection_path(base: Path) -> Path:
    return base / "collection.json"


def get_log_path(base: Path) -> Path:
    return base / "collection.log"


def get_collection_log(venv_path: Path) -> str:
    settings = Inject(get_settings)
    package_path = settings.FRACTAL_ROOT / venv_path  # type: ignore
    log_path = get_log_path(package_path)
    log = log_path.open().read()
    return log


class _TaskCollectPip(TaskCollectPip):
    """
    Internal TaskCollectPip schema

    The difference with its parent class is that we check if the package
    corresponds to a path in the filesystem.
    """

    package_path: Optional[Path] = None

    @property
    def is_local_package(self) -> bool:
        return bool(self.package_path)

    @root_validator(pre=True)
    def check_local_package(cls, values):
        """
        Checks if package corresponds to an existing path on the filesystem

        In this case, the user is providing directly a package file, rather
        than a remote one from PyPI. We set the `package_path` attribute and
        get the actual package name and version from the package file name.
        """
        if "/" in values["package"]:
            package_path = Path(values["package"])
            if not package_path.is_absolute():
                raise ValueError("Package path must be absolute")
            if package_path.exists():
                values["package_path"] = package_path
                (
                    values["package"],
                    values["version"],
                    *_,
                ) = package_path.name.split("-")
        return values

    @property
    def pip_package_version(self):
        """
        Return pip compatible specification of package and version
        """
        version = f"=={self.version}" if self.version else ""
        return f"{self.package}{version}"

    @property
    def source(self):
        if self.is_local_package:
            return f"pip-local:{self.package_path.name}"
        else:
            return f"pip:{self.pip_package_version}"

    def check(self):
        if not self.version:
            raise ValueError("Version is not set or cannot be determined")


def _package_from_path(wheel_path: Path) -> Tuple[str, str]:
    """
    Extract package name and version from package files such as wheel files.
    """
    wheel_filename = wheel_path.name
    package, version, *_rest = wheel_filename.split("-")
    return package, version


def create_package_dir_pip(
    *,
    task_pkg: _TaskCollectPip,
    user: Optional[str] = None,
    **_,
) -> Path:
    settings = Inject(get_settings)
    user = user or settings.FRACTAL_PUBLIC_TASK_SUBDIR

    package_dir = f"{task_pkg.package}{task_pkg.version or ''}"
    venv_path = settings.FRACTAL_ROOT / user / package_dir  # type: ignore
    # TODO check the access right of the venv_path and subdirs
    venv_path.mkdir(exist_ok=False, parents=True)
    return venv_path


async def download_package(
    *,
    task_pkg: _TaskCollectPip,
    dest: Path,
):
    """
    Download package to destination
    """
    interpreter = f"python{task_pkg.python_version}"
    pip = f"{interpreter} -m pip"
    cmd = (
        f"{pip} download --no-deps {task_pkg.pip_package_version} "
        f"-d {dest}"
    )
    stdout = await execute_command(command=cmd, cwd=Path("."))
    pkg_file = next(
        line.split()[-1] for line in stdout.split("\n") if "Saved" in line
    )
    return Path(pkg_file)


def inspect_package(path: Path):
    """
    Inspect task package for version and manifest

    Parameters:
    path: Path
        the path in which the package is saved

    Return
    ------
    version_manifest: dict
        A dictionary containing `version`, the version of the pacakge, and
        `manifest`, the Fractal manifest object relative to the tasks.
    """
    if "whl" in path.as_posix():
        # it is simply a zip file
        # we can extract the version number from *.dist-info/METADATA
        # and read the fractal manifest from the package content
        with ZipFile(path) as wheel:
            namelist = wheel.namelist()
            metadata = next(
                name for name in namelist if "dist-info/METADATA" in name
            )
            manifest = next(
                name for name in namelist if "__FRACTAL_MANIFEST__" in name
            )

            with wheel.open(metadata) as metadata_fd:
                meta = metadata_fd.read().decode("utf-8")
                version = next(
                    line.split()[-1]
                    for line in meta.splitlines()
                    if line.startswith("Version")
                )

            with wheel.open(manifest) as manifest_fd:
                manifest_obj = read_manifest(manifest_fd)  # type: ignore

    version_manifest = dict(version=version, manifest=manifest_obj)
    return version_manifest


async def create_package_environment_pip(
    *,
    task_pkg: _TaskCollectPip,
    venv_path: Path,
) -> List[TaskCreate]:
    """
    Create environment and install package
    """
    logger_name = task_pkg.package
    logger = set_logger(
        logger_name=logger_name,
        log_file_path=get_log_path(venv_path),
        level=logging.DEBUG,
    )
    logger.debug("Creating venv and installing package")

    python_bin, package_root = await _create_venv_install_package(
        path=venv_path,
        task_pkg=task_pkg,
        logger_name=logger_name,
    )

    logger.debug("loading manifest")

    task_list = load_manifest(
        package_root=package_root,
        python_bin=python_bin,
        source=task_pkg.source,
    )
    logger.debug("manifest loaded")
    close_logger(logger)
    return task_list


def read_manifest(file: Union[Path, IOBase]) -> ManifestV1:
    """
    Read and parse manifest file
    """
    if isinstance(file, IOBase):
        manifest_dict = json.load(file)
    else:
        with file.open("r") as f:
            manifest_dict = json.load(f)

    manifest_version = str(manifest_dict["manifest_version"])
    if manifest_version == "1":
        manifest = ManifestV1(**manifest_dict)
    else:
        raise ValueError("Manifest version {manifest_version=} not supported")

    return manifest


def load_manifest(
    package_root: Path,
    python_bin: Path,
    source: str,
) -> List[TaskCreate]:

    manifest_file = package_root / "__FRACTAL_MANIFEST__.json"
    manifest = read_manifest(manifest_file)

    task_list = []
    for t in manifest.task_list:
        task_executable = package_root / t.executable
        if not task_executable.exists():
            raise FileNotFoundError(
                f"Cannot find executable `{task_executable}` "
                f"for task `{t.name}`"
            )
        cmd = f"{python_bin.as_posix()} {task_executable.as_posix()}"
        this_task = TaskCreate(**t.dict(), command=cmd, source=source)
        task_list.append(this_task)
    return task_list


async def _create_venv_install_package(
    *,
    task_pkg: _TaskCollectPip,
    path: Path,
    logger_name: str,
) -> Tuple[Path, Path]:
    """
    Create venv and install package

    Parameters
    ----------
    path : Path
        the directory in which to create the environment
    task_pkg : _TaskCollectPip
        object containing the different metadata required to install the
        package

    Return
    ------
    python_bin: Path
        path to venv's python interpreter
    package_root : Path
        the location of the package manifest
    """
    python_bin = await _init_venv(
        path=path,
        python_version=task_pkg.python_version,
        logger_name=logger_name,
    )
    package_root = await _pip_install(
        venv_path=path, task_pkg=task_pkg, logger_name=logger_name
    )
    return python_bin, package_root


async def _init_venv(
    *,
    path: Path,
    python_version: str = "3.8",
    logger_name: str,
) -> Path:
    """
    Set a virtual environment at `path/venv`

    Parameters
    ----------
    path : Path
        path to directory in which to set up the virtual environment
    python_version : str, default='3.8'
        Python version the virtual environment will be based upon

    Return
    ------
    python_bin : Path
        path to python interpreter
    """
    interpreter = f"python{python_version}"
    await execute_command(
        cwd=path, command=f"{interpreter} -m venv venv", logger_name="fractal"
    )
    return path / "venv/bin/python"


async def _pip_install(
    venv_path: Path,
    task_pkg: _TaskCollectPip,
    logger_name: str,
) -> Path:
    """
    Install package in venv

    Return
    ------
    package_root : Path
        the location of the package manifest
    """
    pip = venv_path / "venv/bin/pip"

    if task_pkg.is_local_package:
        pip_install_str = task_pkg.package_path.as_posix()  # type: ignore
    else:
        version_string = f"=={task_pkg.version}" if task_pkg.version else ""
        extras = (
            f"[{task_pkg.package_extras}]" if task_pkg.package_extras else ""
        )
        pip_install_str = f"{task_pkg.package}{extras}{version_string}"

    cmd_install = f"{pip} install {pip_install_str}"
    cmd_inspect = f"{pip} show {task_pkg.package}"

    await execute_command(
        cwd=venv_path,
        command=f"{pip} install --upgrade pip",
        logger_name=logger_name,
    )
    await execute_command(
        cwd=venv_path, command=cmd_install, logger_name=logger_name
    )

    # Extract package installation path from `pip show`
    stdout_inspect = await execute_command(
        cwd=venv_path, command=cmd_inspect, logger_name=logger_name
    )

    location = Path(
        next(
            line.split()[-1]
            for line in stdout_inspect.split("\n")
            if line.startswith("Location:")
        )
    )
    package_root = location / task_pkg.package.replace("-", "_")
    if not package_root.exists():
        raise RuntimeError(
            "Could not determine package installation location."
        )
    return package_root
