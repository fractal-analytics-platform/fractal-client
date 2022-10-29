"""
Copyright 2022 (C) Friedrich Miescher Institute for Biomedical Research and
University of Zurich

Original authors:
Jacopo Nespolo <jacopo.nespolo@exact-lab.it>

This file is part of Fractal and was originally developed by eXact lab S.r.l.
<exact-lab.it> under contract with Liberali Lab from the Friedrich Miescher
Institute for Biomedical Research and Pelkmans Lab from the University of
Zurich.
"""
import json
from pathlib import Path
from typing import List
from typing import Literal
from typing import Optional
from typing import Tuple

from pydantic import root_validator

from ..app.schemas import ManifestV1
from ..app.schemas import TaskCollectPip
from ..app.schemas import TaskCreate
from ..config import get_settings
from ..syringe import Inject
from ..utils import execute_command


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
        package_path = Path(values["package"])
        if package_path.exists():
            values["package_path"] = package_path
            values["package"], values["version"], *_ = package_path.name.split(
                "-"
            )
        return values

    @property
    def source(self):
        if self.is_local_package:
            return f"pip-local:{self.package_path.name}"
        else:
            return f"pip:{self.package}=={self.version}"


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
    user: str = ".fractal",
    **_,
) -> Path:
    settings = Inject(get_settings)
    # assemble installation path
    package_path = Path(task_pkg.package)
    if package_path.exists():
        # The package is a local package in the filesystem rather than a
        # remote PyPI package. Extract package and version from basename
        package, version = _package_from_path(package_path)

    package_dir = f"{task_pkg.package}{task_pkg.version or ''}"
    if settings.FRACTAL_ROOT:
        # It should be always true, we are only checking that the system is
        # fully configured
        env_path = settings.FRACTAL_ROOT / user / package_dir
    # TODO check the access right of the env_path and subdirs
    env_path.mkdir(exist_ok=False, parents=True)
    return env_path


async def create_package_environment_pip(
    *,
    task_pkg: _TaskCollectPip,
    env_path: Path,
    user: str = ".fractal",
    env_type: Literal["venv"] = "venv",
) -> List[TaskCreate]:
    """
    Create environment and install package
    """

    if env_type == "venv":
        python_bin, package_root = await _create_venv_install_package(
            path=env_path,
            task_pkg=task_pkg,
        )
    else:
        raise ValueError(f"Environment type {env_type} not supported")

    task_list = load_manifest(
        package_root=package_root,
        python_bin=python_bin,
        source=task_pkg.source,
    )
    return task_list


def load_manifest(
    package_root: Path,
    python_bin: Path,
    source: str,
) -> List[TaskCreate]:

    manifest_file = package_root / "__FRACTAL_MANIFEST__.json"
    with manifest_file.open("r") as f:
        manifest_dict = json.load(f)

    task_list = []
    if str(manifest_dict["manifest_version"]) == "1":
        manifest = ManifestV1(**manifest_dict)

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
        path=path, python_version=task_pkg.python_version
    )
    package_root = await _pip_install(venv_path=path, task_pkg=task_pkg)
    return python_bin, package_root


async def _init_venv(*, path: Path, python_version: str = "3.8") -> Path:
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
    cmd_inspect = f"{pip} show -f {task_pkg.package}"

    await execute_command(
        cwd=venv_path, command=cmd_install, logger_name="fractal"
    )
    stdout_inspect = await execute_command(
        cwd=venv_path, command=cmd_inspect, logger_name="fractal"
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
