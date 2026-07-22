import argparse as ap


def add_task_parser(subparsers):
    task_parser = subparsers.add_parser(
        "task",
        description="Task commands.",
        allow_abbrev=False,
    )
    task_subparsers = task_parser.add_subparsers(
        title="Valid sub-commands", dest="subcmd", required=True
    )

    # task list
    task_list_parser = task_subparsers.add_parser(
        "list",
        description="List tasks.",
        allow_abbrev=False,
    )

    # task make-core
    task_make_core_parser = task_subparsers.add_parser(
        "make-core",
        description="Mark tasks as core.",
        allow_abbrev=False,
    )
    task_make_core_parser.add_argument(
        "task_ids",
        help="List of the IDs of tasks that should be made core.",
        type=int,
        nargs="+",
    )

    # task make-not-core
    task_make_not_core_parser = task_subparsers.add_parser(
        "make-not-core",
        description="Mark tasks as non-core.",
        allow_abbrev=False,
    )
    task_make_not_core_parser.add_argument(
        "task_ids",
        help="List of the IDs of tasks that should be made core.",
        type=int,
        nargs="+",
    )
    # task collect
    task_collect_parser = task_subparsers.add_parser(
        "collect",
        description=(
            "Install and collect all tasks from a pip-installable package."
        ),
        allow_abbrev=False,
    )
    task_collect_parser.add_argument(
        "package",
        help="Package name or absolute path to local package.",
    )
    task_collect_parser.add_argument(
        "--python-version",
        help="Select the python version to use for this package.",
    )
    task_collect_parser.add_argument(
        "--package-version",
        help="Select the package version.",
    )
    task_collect_parser.add_argument(
        "--package-extras",
        help=(
            "Comma separated list of extra components for the package to be "
            "installed, e.g., `collect fractal-tasks-core "
            "--package-extras=torch,tensorflow` will trigger the installation of "
            "`fractal-tasks-core[torch,tensorflow]`."
        ),
    )
    task_collect_parser.add_argument(
        "--pre-pinned-dependency",
        action="append",
        help=(
            "Package/version pair representing a pre-pinned-version dependency, "
            "in the form `collect fractal-tasks-core --pre-pinned-dependency "
            "pydantic=1.10.0`. Include `--pre-pinned-dependency` multiple times "
            "to pin several packages to specific versions."
        ),
    )
    task_collect_parser.add_argument(
        "--post-pinned-dependency",
        action="append",
        help=(
            "Package/version pair representing a post-pinned-version dependency, "
            "in the form `collect fractal-tasks-core --post-pinned-dependency "
            "pydantic=1.10.0`. Include `--post-pinned-dependency` multiple times "
            "to pin several packages to specific versions."
        ),
    )
    task_collect_parser.add_argument(
        "--private",
        default=False,
        action="store_true",
        help="Make task group private.",
    )

    # task collect custom
    task_collect_custom_parser = task_subparsers.add_parser(
        "collect-custom",
        description="Collect all tasks from a custom Python interpreter.",
        allow_abbrev=False,
    )
    task_collect_custom_parser.add_argument(
        "label",
        help="A common label identifying this package.",
    )
    task_collect_custom_parser.add_argument(
        "python_interpreter",
        help=(
            "Absolute path to the Python interpreter to be used for running tasks."
        ),
    )
    task_collect_custom_parser.add_argument(
        "manifest",
        help="Local path of the Manifest of the Fractal task package.",
    )
    task_collect_custom_parser.add_argument(
        "--version",
        help="Version of tasks to be collected.",
        required=True,
    )
    tasktask_collect_custom_pkg_name_or_root = (
        task_collect_custom_parser.add_mutually_exclusive_group(required=True)
    )
    tasktask_collect_custom_pkg_name_or_root.add_argument(
        "--package-name",
        help=(
            "Name of the package, as used in 'import <package_name>'; "
            "this is then used to extract the package directory (package_root) "
            "via 'importlib.util.find_spec <package_name>'."
        ),
    )
    tasktask_collect_custom_pkg_name_or_root.add_argument(
        "--package-root",
        help=(
            "The folder where the package is installed. If not provided, "
            "it will be  automatically inferred based on 'package_name'."
        ),
    )
    task_collect_custom_parser.add_argument(
        "--private",
        default=False,
        action="store_true",
        help="Make task group private.",
    )

    # task check-collection
    task_check_collection_parser = task_subparsers.add_parser(
        "check-collection",
        description="Check status of background task collection processes.",
        allow_abbrev=False,
    )
    task_check_collection_parser.add_argument(
        "task_group_activity_id",
        help="Activity ID of the collection (see output of `task collect`).",
        type=int,
    )
    task_check_collection_parser.add_argument(
        "--include-logs",
        default=False,
        action="store_true",
        help="Also include task-collection logs.",
    )

    # task new
    task_new_parser = task_subparsers.add_parser(
        "new",
        description="Create new task.",
        argument_default=ap.SUPPRESS,
        allow_abbrev=False,
    )
    task_new_parser.add_argument(
        "name", help="A human readable name for the task."
    )
    task_new_parser.add_argument(
        "--task-type",
        help="The task type (e.g. 'parallel', 'non_parallel', 'compound').",
    )
    task_new_parser.add_argument(
        "--command-non-parallel",
        help="The non parallel command that executes the task.",
    )
    task_new_parser.add_argument(
        "--command-parallel",
        help="The  parallel command that executes the task.",
    )
    task_new_parser.add_argument(
        "--version",
        help="Task version.",
        required=True,
    )
    task_new_parser.add_argument(
        "--meta-non-parallel",
        help="Path to JSON file with meta non parallel arguments.",
    )
    task_new_parser.add_argument(
        "--meta-parallel",
        help="Path to JSON file with meta parallel arguments.",
    )
    task_new_parser.add_argument(
        "--args-schema-non-parallel",
        help="Path to JSON file with args non parallel arguments.",
    )
    task_new_parser.add_argument(
        "--args-schema-parallel",
        help="Path to JSON file with arg parallel arguments.",
    )
    task_new_parser.add_argument(
        "--args-schema-version",
        help=(
            "Label encoding how the task-arguments JSON Schema was generated "
            "(e.g. `pydantic_v1`)."
        ),
    )
    task_new_parser.add_argument(
        "--private",
        default=False,
        action="store_true",
        help="Make task group private.",
    )

    # task edit
    task_edit_parser = task_subparsers.add_parser(
        "edit",
        description="Edit task.",
        argument_default=ap.SUPPRESS,
        allow_abbrev=False,
    )

    task_edit_id_or_name_group = task_edit_parser.add_mutually_exclusive_group(
        required=True
    )
    task_edit_id_or_name_group.add_argument(
        "--id", help="ID of the task to edit.", type=int
    )
    task_edit_id_or_name_group.add_argument(
        "--name", help="Name of the task to edit."
    )

    task_edit_parser.add_argument(
        "--version",
        help=(
            "Version of the task to edit "
            "(only accepted in combination with `--name`)."
        ),
    )
    task_edit_parser.add_argument(
        "--command-non-parallel", help="New task non parallel command."
    )
    task_edit_parser.add_argument(
        "--command-parallel",
        help="New task parallel command.",
    )
    task_edit_parser.add_argument(
        "--input-types",
        help=("Path to JSON file with new input types."),
    )
    task_edit_parser.add_argument(
        "--output-types",
        help=("Path to JSON file with new output types."),
    )
