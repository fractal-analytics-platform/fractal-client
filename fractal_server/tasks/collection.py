import asyncio
import json
from pathlib import Path
from shlex import split as shlex_split
from typing import List
from typing import Literal
from typing import Optional
from typing import Tuple

from ..app.schemas import ManifestV1
from ..app.schemas import TaskCreate
from ..config import get_settings
from ..syringe import Inject


async def _execute_command(*, cwd: Path, command: str) -> str:
    command_split = shlex_split(command)
    cmd, *args = command_split

    proc = await asyncio.create_subprocess_exec(
        cmd,
        *args,
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.PIPE,
        cwd=cwd,
    )
    stdout, stderr = await proc.communicate()
    if proc.returncode != 0:
        raise RuntimeError(stderr.decode("utf-8"))
    return stdout.decode("utf-8")


async def create_package_environment(
    *,
    package: str,
    version: Optional[str],
    user: str = ".fractal",
    python_version: str = "3.8",
    env_type: Literal["venv"] = "venv",
) -> List[TaskCreate]:
    """
    Create environment and install package
    """
    settings = Inject(get_settings)
    # assemble installation path
    package_dir = f"{package}{version or ''}"
    if settings.FRACTAL_ROOT:
        env_path = settings.FRACTAL_ROOT / user / package_dir
    # TODO check the access right of the env_path and subdirs
    env_path.mkdir(exist_ok=True, parents=True)

    if env_type == "venv":
        python_bin, package_root = await _create_venv_install_package(
            path=env_path,
            package=package,
            version=version,
            python_version=python_version,
        )
    else:
        raise ValueError(f"Environment type {env_type} not supported")

    source = f"pypi:{package}=={version}"
    task_list = load_manifest(
        package_root=package_root,
        python_bin=python_bin,
        source=source,
    )
    return task_list


def load_manifest(
    package_root: Path,
    python_bin: Path,
    source: str,
) -> List[TaskCreate]:

    manifest_file = package_root / "__FRACTAL__MANIFEST__.json"
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
    path: Path,
    package: str,
    version: Optional[str],
    python_version: str,
) -> Tuple[Path, Path]:
    """
    Create venv and install package

    Parameters
    ----------
    path : Path
        the directory in which to create the environment
    package : str
        package name
    version : str, optional
        package version
    python_version : str
        version of the Python interpreter

    Return
    ------
    python_bin: Path
        path to venv's python interpreter
    package_root : Path
        the location of the package manifest
    """
    python_bin = await _init_venv(path=path, python_version=python_version)
    package_root = await _pip_install(
        venv_path=path, package=package, version=version
    )
    return python_bin, package_root


async def _init_venv(*, path: Path, python_version: str) -> Path:
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
    await _execute_command(cwd=path, command=f"{interpreter} -m venv venv")
    return path / "venv/bin/python"


async def _pip_install(
    venv_path: Path, package: str, version: Optional[str]
) -> Path:
    """
    Install package in venv

    Return
    ------
    package_root : Path
        the location of the package manifest
    """
    pip = venv_path / "venv/bin/pip"
    version_string = f"=={version}" if version else ""

    cmd_install = f"{pip} install {package}{version_string}"
    cmd_inspect = f"{pip} show -f {package}"

    await _execute_command(cwd=venv_path, command=cmd_install)
    stdout_inspect = await _execute_command(cwd=venv_path, command=cmd_inspect)

    location = Path(
        next(
            line.split()[-1]
            for line in stdout_inspect.split("\n")
            if line.startswith("Location:")
        )
    )
    package_root = location / package.replace("-", "_")
    if not package_root.exists():
        raise RuntimeError(
            "Could not determine package installation location."
        )
    return package_root
