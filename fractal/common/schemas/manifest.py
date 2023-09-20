from typing import Any
from typing import Optional
from typing import TypeVar

from pydantic import BaseModel
from pydantic import Field
from pydantic import HttpUrl
from pydantic import root_validator
from pydantic import validator


__all__ = ("TaskManifestV1", "ManifestV1")


class _TaskManifestBase(BaseModel):
    """
    Base class for `TaskManifest`

    Represents a task within a manfest

    Attributes:
        name:
            The task name
        executable:
            Path to the executable relative to the package root

            Note: by package root we mean "as it will be installed". If a
            package `Pkg` installs in the folder `pkg` the executable
            `pkg/executable.py`, this attribute must contain only
            `executable.py`.
        input_type:
            The input type accepted by the task
        output_type:
            The output type returned by the task
        meta:
            Additional information about the package, such as hash of the
            executable, specific runtime requirements (e.g., need_gpu=True),
            etc.
        args_schema:
            JSON Schema for task arguments
        docs_info:
            Additional information about the Task, coming from the docstring.
        docs_link:
            Link to Task docs.
    """

    name: str
    executable: str
    input_type: str
    output_type: str
    meta: Optional[dict[str, Any]] = Field(default_factory=dict)
    args_schema: Optional[dict[str, Any]]
    docs_info: Optional[str]
    docs_link: Optional[HttpUrl]


TaskManifestType = TypeVar("TaskManifestType", bound=_TaskManifestBase)


class _ManifestBase(BaseModel):
    """
    Manifest base class

    Packages containing tasks are required to include a special file
    `__FRACTAL_MANIFEST__.json` in order to be discovered and used by Fractal.

    This model class and the model classes it depends on provide the base
    schema to read, write and validate manifests.

    Attributes:
        manifest_version:
            A version string that provides indication for compatibility between
            manifests as the schema evolves. This is for instance used by
            Fractal to determine which subclass of the present base class needs
            be used to read and validate the input.
        task_list : list[TaskManifestType]
            The list of tasks, represented as specified by subclasses of the
            _TaskManifestBase (a.k.a. TaskManifestType)
        has_args_schemas:
            `True` if the manifest incldues JSON Schemas for the arguments of
            each task.
        args_schema_version:
            Label of how `args_schema`s were generated (e.g. `pydantic_v1`).
    """

    manifest_version: str
    task_list: list[TaskManifestType]
    has_args_schemas: bool = False
    args_schema_version: Optional[str]

    @root_validator()
    def _check_args_schemas_are_present(cls, values):
        has_args_schemas = values["has_args_schemas"]
        task_list = values["task_list"]
        if has_args_schemas:
            for task in task_list:
                if task.args_schema is None:
                    raise ValueError(
                        f'has_args_schemas={has_args_schemas} but task "'
                        f'{task.name}" has args_schema={task.args_schema}.'
                    )
        return values


class TaskManifestV1(_TaskManifestBase):
    pass


class ManifestV1(_ManifestBase):
    """
    Manifest schema version 1
    """

    task_list: list[TaskManifestV1]

    @validator("manifest_version")
    def manifest_version_1(cls, value):
        if value != "1":
            raise ValueError(f"Wrong manifest version (given {value})")
