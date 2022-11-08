from pathlib import Path
from typing import Any
from typing import Dict
from typing import List
from typing import Optional
from typing import TypeVar

from pydantic import BaseModel
from pydantic import validator
from sqlmodel import Field


__all__ = ("TaskManifestV1", "ManifestV1")


class _TaskManifestBase(BaseModel):
    """
    Base class for TaskManifest

    Represents a task within a manfest

    Attributes
    ----------
    name : str
        the task name
    executable : Path
        path to the executable relative to the package root

        Note: by package root we mean "as it will be installed". If a package
        `Pkg` installs in the folder `pkg` the executable `pkg/executable.py`,
        this attribute must contain only `executable.py`.
    input_type : str
        the input type accepted by the task
    output_type : str
        the output type returned by the task
    default_args : Dict[str, Any]
        an arbitrary, JSON-serializable dictionary containing the default
        parameters that will be passed to the task
    meta : Dict[str, Any]
        additional information about the package, such as hash of the
        executable, specific runtime requirements (e.g., need_gpu=True), etc.
    """

    name: str
    executable: Path
    input_type: str
    output_type: str
    default_args: Optional[Dict[str, Any]] = Field(default_factory=dict)
    meta: Optional[Dict[str, Any]] = Field(default_factory=dict)


TaskManifestType = TypeVar("TaskManifestType", bound=_TaskManifestBase)


class _ManifestBase(BaseModel):
    """
    Manifest base class

    Packages containing tasks are required to include a special file
    `__FRACTAL_MANIFEST__.json` in order to be discovered and used by Fractal.

    This model class and the model classes it depends on provide the base
    schema to read, write and validate manifests.

    Attributes
    ----------
    manifest_version : str
        a version string that provides indication for compatibility between
        manifests as the schema evolves. This is for instance used by Fractal
        to determine which subclass of the present base class needs be used to
        read and validate the input.
    task_list : List[TaskManifestType]
        the list of tasks, represented as specified by subclasses of the
        _TaskManifestBase (a.k.a. TaskManifestType)
    """

    manifest_version: str
    task_list: List[TaskManifestType]  # type: ignore


class TaskManifestV1(_TaskManifestBase):
    pass


class ManifestV1(_ManifestBase):
    """
    Manifest schema version 1
    """

    task_list: List[TaskManifestV1]

    @validator("manifest_version")
    def manifest_version_1(cls, value):
        if value != "1":
            raise ValueError("Wrong manifest version")
