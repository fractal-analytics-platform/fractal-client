"""
Copyright 2022 (C) Friedrich Miescher Institute for Biomedical Research and
University of Zurich

Original authors:
Jacopo Nespolo <jacopo.nespolo@exact-lab.it>
Marco Franzon <marco.franzon@exact-lab.it>
Tommaso Comparin <tommaso.comparin@exact-lab.it>
Yuri Chiucconi  <yuri.chiucconi@exact-lab.it>

This file is part of Fractal and was originally developed by eXact lab S.r.l.
<exact-lab.it> under contract with Liberali Lab from the Friedrich Miescher
Institute for Biomedical Research and Pelkmans Lab from the University of
Zurich.
"""
import argparse as ap


parser_main = ap.ArgumentParser(
    description="Command-line interface for Fractal Client.",
    allow_abbrev=False,
)

parser_main.add_argument(
    "-u",
    "--user",
    help=(
        "User email address for login (overrides `FRACTAL_USER` "
        "environment variable)."
    ),
)
parser_main.add_argument(
    "-p",
    "--password",
    help="User password (overrides `FRACTAL_PASSWORD` environment variable).",
)

parser_main.add_argument(
    "--token-path",
    help="User token (overrides `FRACTAL_TOKEN_PATH` environment variable).",
)
parser_main.add_argument(
    "--fractal-server",
    help=(
        "URL of Fractal server (overrides `FRACTAL_SERVER` "
        "environment variable)."
    ),
)

parser_main.add_argument(
    "--debug",
    action="store_true",
    default=False,
    help="Change minimal logging level from INFO to DEBUG.",
)
parser_main.add_argument(
    "--batch",
    default=False,
    action="store_true",
    help=(
        "Return output suitable for scripting, e.g., "
        "only the ID of items created instead of the full object."
    ),
)

subparsers_main = parser_main.add_subparsers(title="Commands", dest="cmd")


# PROJECT GROUP
project_parser = subparsers_main.add_parser(
    "project",
    description="Project commands.",
    allow_abbrev=False,
)
project_subparsers = project_parser.add_subparsers(
    title="Valid sub-commands", dest="subcmd", required=True
)

# project new
project_new_parser = project_subparsers.add_parser(
    "new",
    description="Create new project.",
    allow_abbrev=False,
)
project_new_parser.add_argument("name", help="Name of new project.")

# project list
project_list_parser = project_subparsers.add_parser(
    "list",
    description="List projects.",
    allow_abbrev=False,
)

# project show
project_show_parser = project_subparsers.add_parser(
    "show",
    description="Show details of single project.",
    allow_abbrev=False,
)
project_show_parser.add_argument(
    "project_id", type=int, help="ID of project to show."
)

# project delete
project_delete_parser = project_subparsers.add_parser(
    "delete",
    description="Delete project.",
    allow_abbrev=False,
)
project_delete_parser.add_argument(
    "project_id", type=int, help="ID of project to delete."
)

# project add-dataset
project_add_dataset_parser = project_subparsers.add_parser(
    "add-dataset",
    description="Add dataset to project.",
    allow_abbrev=False,
)
project_add_dataset_parser.add_argument(
    "project_id",
    type=int,
    help="ID of project to add the new dataset to.",
)
project_add_dataset_parser.add_argument(
    "dataset_name",
    help="Name of new dataset.",
)
project_add_dataset_parser.add_argument(
    "--zarr-dir",
    help="Path to zarr dir.",
    required=False,
)


# project edit
project_edit_parser = project_subparsers.add_parser(
    "edit",
    description="Edit details of a single project.",
    allow_abbrev=False,
)
project_edit_parser.add_argument(
    "project_id", type=int, help="ID of the project to edit."
)
project_edit_parser.add_argument(
    "--new-name", help="New project name.", required=False
)


# DATASET GROUP
dataset_parser = subparsers_main.add_parser(
    "dataset",
    description="Dataset commands.",
    allow_abbrev=False,
)
dataset_subparsers = dataset_parser.add_subparsers(
    title="Valid sub-commands", dest="subcmd", required=True
)


# dataset edit
dataset_edit_parser = dataset_subparsers.add_parser(
    "edit",
    description="Edit dataset.",
    argument_default=ap.SUPPRESS,
    allow_abbrev=False,
)
dataset_edit_parser.add_argument("project_id", type=int, help="Project ID.")
dataset_edit_parser.add_argument("dataset_id", type=int, help="Dataset ID.")
dataset_edit_parser.add_argument("--new-name", help="New name of dataset.")


# dataset show
dataset_show_parser = dataset_subparsers.add_parser(
    "show",
    description="Show dataset.",
    argument_default=ap.SUPPRESS,
    allow_abbrev=False,
)
dataset_show_parser.add_argument("project_id", type=int, help="Project ID.")
dataset_show_parser.add_argument("dataset_id", type=int, help="Dataset ID.")

# dataset delete
dataset_delete_parser = dataset_subparsers.add_parser(
    "delete",
    description="Delete dataset.",
    argument_default=ap.SUPPRESS,
    allow_abbrev=False,
)
dataset_delete_parser.add_argument("project_id", type=int, help="Project ID.")
dataset_delete_parser.add_argument("dataset_id", type=int, help="Dataset ID.")

# TASK GROUP
task_parser = subparsers_main.add_parser(
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
    "--pinned-dependency",
    action="append",
    help=(
        "Package/version pair representing a pinned-version dependency, in "
        "the form `collect fractal-tasks-core --pinned-dependency "
        "pydantic=1.10.0`. Include `--pinned-dependency` multiple times to "
        "pin several packages to specific versions."
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
    "manifest", help="Local path of the Manifest of the Fractal task package."
)
task_collect_custom_parser.add_argument(
    "--version",
    help="Version of tasks to be collected.",
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
    "--command-parallel", help="The  parallel command that executes the task."
)
task_new_parser.add_argument(
    "--version",
    help="Task version.",
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
task_edit_parser.add_argument("--new-version", help="New task version.")
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

# WORKFLOW GROUP

workflow_parser = subparsers_main.add_parser(
    "workflow",
    description="Workflow commands.",
    allow_abbrev=False,
)
workflow_subparsers = workflow_parser.add_subparsers(
    title="Valid sub-commands", dest="subcmd", required=True
)

# workflow list
workflow_list_parser = workflow_subparsers.add_parser(
    "list",
    description="List workflows for given project.",
    allow_abbrev=False,
)
workflow_list_parser.add_argument("project_id", type=int, help="Project ID.")

# workflow new
workflow_new_parser = workflow_subparsers.add_parser(
    "new",
    description="Create new workflow.",
    allow_abbrev=False,
)
workflow_new_parser.add_argument(
    "name",
    help="Workflow name (must be unique, and not only made of numbers only).",
)
workflow_new_parser.add_argument("project_id", type=int, help="Project ID.")

# workflow show
workflow_show_parser = workflow_subparsers.add_parser(
    "show",
    description="Show workflow.",
    allow_abbrev=False,
)
workflow_show_parser.add_argument("project_id", type=int, help="Project ID.")
workflow_show_parser.add_argument("workflow_id", type=int, help="Workflow ID.")

# workflow edit
workflow_edit_parser = workflow_subparsers.add_parser(
    "edit",
    description="Edit workflow.",
    argument_default=ap.SUPPRESS,
    allow_abbrev=False,
)
workflow_edit_parser.add_argument("project_id", type=int, help="Project ID.")
workflow_edit_parser.add_argument("workflow_id", type=int, help="Workflow ID.")
workflow_edit_parser.add_argument(
    "--new-name",
    help="New workflow name.",
    required=True,
)

# workflow delete
workflow_delete_parser = workflow_subparsers.add_parser(
    "delete",
    description="Delete workflow.",
    allow_abbrev=False,
)
workflow_delete_parser.add_argument("project_id", type=int, help="Project ID.")
workflow_delete_parser.add_argument(
    "workflow_id", type=int, help="Workflow ID."
)


# workflow add task
workflow_add_task_parser = workflow_subparsers.add_parser(
    "add-task",
    description="Add new task to specific workflow.",
    allow_abbrev=False,
)
workflow_add_task_parser.add_argument(
    "project_id", type=int, help="Project ID."
)
workflow_add_task_parser.add_argument(
    "workflow_id",
    type=int,
    help="Workflow ID.",
)

workflow_add_task_id_or_name_group = (
    workflow_add_task_parser.add_mutually_exclusive_group(required=True)
)
workflow_add_task_id_or_name_group.add_argument(
    "--task-id", help="ID of the task to add.", type=int
)
workflow_add_task_id_or_name_group.add_argument(
    "--task-name", help="Name of the task to add."
)
workflow_add_task_parser.add_argument(
    "--task-version",
    help=(
        "Version of task to add "
        "(only accepted in combination with --task-name)."
    ),
)
workflow_add_task_parser.add_argument(
    "--order", help="Order of this task within the workflow's task list."
)

workflow_add_task_parser.add_argument(
    "--args-non-parallel", help="Args for non parallel tasks"
)

workflow_add_task_parser.add_argument(
    "--args-parallel", help="Args for parallel tasks"
)

workflow_add_task_parser.add_argument(
    "--meta-non-parallel", help="Metadata file for non-parallel tasks"
)

workflow_add_task_parser.add_argument(
    "--meta-parallel", help="Metadata file for parallel tasks"
)

workflow_add_task_parser.add_argument(
    "--type-filters",
    help="Path to JSON file with type filters.",
)


# workflow edit task
workflow_edit_task_parser = workflow_subparsers.add_parser(
    "edit-task",
    description="Edit task within specific workflow.",
    allow_abbrev=False,
)
workflow_edit_task_parser.add_argument(
    "project_id", type=int, help="Project ID."
)
workflow_edit_task_parser.add_argument(
    "workflow_id",
    type=int,
    help="Workflow ID.",
)
workflow_edit_task_parser.add_argument(
    "workflow_task_id",
    type=int,
    help="Workflow task ID, the ID of a task inside the list of tasks.",
)
workflow_edit_task_parser.add_argument(
    "--type-filters",
    help=(
        "Path to JSON file containing the type filters of this "
        "workflow task."
    ),
)
workflow_edit_task_parser.add_argument(
    "--args-non-parallel", help="Args for non parallel tasks"
)

workflow_edit_task_parser.add_argument(
    "--args-parallel", help="Args for parallel tasks"
)

workflow_edit_task_parser.add_argument(
    "--meta-non-parallel", help="Metadata file fornon parallel tasks"
)

workflow_edit_task_parser.add_argument(
    "--meta-parallel", help="Metadata file for parallel tasks"
)


# workflow remove task
workflow_remove_task_parser = workflow_subparsers.add_parser(
    "rm-task",
    description="Remove task from a specific workflow.",
    allow_abbrev=False,
)
workflow_remove_task_parser.add_argument(
    "project_id", type=int, help="Project ID."
)
workflow_remove_task_parser.add_argument(
    "workflow_id",
    type=int,
    help="Workflow ID.",
)
workflow_remove_task_parser.add_argument(
    "workflow_task_id",
    type=int,
    help="Workflow task ID (the ID of a task inside the list of tasks).",
)


# workflow import
workflow_import_parser = workflow_subparsers.add_parser(
    "import",
    description="Import workflow to project from file.",
    allow_abbrev=False,
)
workflow_import_parser.add_argument(
    "--project-id",
    type=int,
    help="ID of the project where the workflow will be imported.",
    required=True,
)
workflow_import_parser.add_argument(
    "--json-file",
    type=str,
    help="Path to a JSON file with the workflow to be imported.",
    required=True,
)
workflow_import_parser.add_argument(
    "--workflow-name",
    type=str,
    help="Name of the new workflow (if set, overrides the one in JSON file)",
    required=False,
)

# workflow export
workflow_export_parser = workflow_subparsers.add_parser(
    "export",
    description="Export workflow to file.",
    allow_abbrev=False,
)
workflow_export_parser.add_argument(
    "project_id",
    type=int,
    help="Project ID.",
)
workflow_export_parser.add_argument(
    "workflow_id",
    type=int,
    help="Workflow ID.",
)
workflow_export_parser.add_argument(
    "--json-file",
    help="Path to the JSON file where the workflow will be exported.",
    required=True,
)

# JOB GROUP

job_parser = subparsers_main.add_parser(
    "job",
    description="Job commands.",
    allow_abbrev=False,
)
job_subparsers = job_parser.add_subparsers(
    title="Valid sub-commands", dest="subcmd", required=True
)

# job list
job_list_parser = job_subparsers.add_parser(
    "list",
    description="List jobs for given project.",
    allow_abbrev=False,
)
job_list_parser.add_argument(
    "project_id",
    type=int,
    help="Project ID.",
)

# job show
job_show_parser = job_subparsers.add_parser(
    "show",
    description="Query status of workflow-execution job.",
    argument_default=ap.SUPPRESS,
    allow_abbrev=False,
)
job_show_parser.add_argument("project_id", type=int, help="Project ID.")
job_show_parser.add_argument("job_id", type=int, help="Job ID.")

# job download-logs
job_download_logs_parser = job_subparsers.add_parser(
    "download-logs",
    description="Download full folder of workflow-execution job.",
    allow_abbrev=False,
)
job_download_logs_parser.add_argument(
    "project_id", type=int, help="Project ID."
)
job_download_logs_parser.add_argument("job_id", type=int, help="Job ID.")
job_download_logs_parser.add_argument(
    "--output",
    dest="output_folder",
    help="Path of the output folder.",
    required=True,
)

# job stop
job_stop_parser = job_subparsers.add_parser(
    "stop",
    description="Stop workflow-execution job.",
    allow_abbrev=False,
)
job_stop_parser.add_argument("project_id", type=int, help="Project ID.")
job_stop_parser.add_argument("job_id", type=int, help="Job ID.")


# job submit
job_submit_parser = job_subparsers.add_parser(
    "submit",
    description="Submit a job.",
    argument_default=ap.SUPPRESS,
    allow_abbrev=False,
)

job_submit_parser.add_argument("project_id", type=int)
job_submit_parser.add_argument("workflow_id", type=int)
job_submit_parser.add_argument("dataset_id", type=int)
job_submit_parser.add_argument(
    "--start",
    dest="first_task_index",
    type=int,
    help=(
        "Positional index of the first task to be executed"
        " (starting from 0)."
    ),
    required=False,
)
job_submit_parser.add_argument(
    "--end",
    dest="last_task_index",
    type=int,
    help=(
        "Positional index of the last task to be executed"
        " (starting from 0)."
    ),
    required=False,
)
job_submit_parser.add_argument(
    "-w",
    "--worker-init",
    help="Command to be run before starting a worker.",
)
job_submit_parser.add_argument(
    "--attribute-filters-json",
    help=(
        "Path to JSON file with the attribute filters for this job submission."
    ),
    required=False,
)
job_submit_parser.add_argument(
    "--type-filters-json",
    help=("Path to JSON file with the type filters for this job submission."),
    required=False,
)


# VERSION GROUP
version_parser = subparsers_main.add_parser(
    "version",
    description="Print version and exit.",
    allow_abbrev=False,
)


# USER GROUP

user_parser = subparsers_main.add_parser(
    "user",
    description="User commands.",
    allow_abbrev=False,
)
user_subparsers = user_parser.add_subparsers(
    title="Valid sub-commands", dest="subcmd", required=True
)

# user whoami
user_whoami_parser = user_subparsers.add_parser(
    "whoami",
    description="Get info on current user (fails if user is not registered).",
    allow_abbrev=False,
)
user_whoami_parser.add_argument(
    "--viewer-paths",
    help="Include user's `viewer_paths` attribute.",
    action="store_true",
    required=False,
)

# user register
user_register_parser = user_subparsers.add_parser(
    "register",
    description=(
        "Register a new user with the Fractal server and edit their settings "
        "(note: user creation and settings editing are two independent steps)."
    ),
    allow_abbrev=False,
)
user_register_parser.add_argument(
    "new_email", help="Email to be used as username."
)
user_register_parser.add_argument(
    "new_password", help="Password for the new user."
)
user_register_parser.add_argument(
    "--project-dir",
    help="User-writeable base folder, used e.g. for default `zarr_dir` paths.",
    required=False,
)
user_register_parser.add_argument(
    "--slurm-user",
    help="Username to login into SLURM cluster.",
    required=False,
)
user_register_parser.add_argument(
    "--username",
    help="Username associated to the user.",
    required=False,
)
user_register_parser.add_argument(
    "--ssh-settings-json",
    help=(
        "Path to JSON file with (a subset of) following settings: "
        "ssh_host, ssh_username, ssh_private_key_path, "
        "ssh_tasks_dir, ssh_jobs_dir."
    ),
    type=str,
    required=False,
)
user_register_parser.add_argument(
    "--superuser",
    help="Give superuser privileges to the new user.",
    action="store_true",
    required=False,
)

# user list
user_list_parser = user_subparsers.add_parser(
    "list",
    description="List all users.",
    allow_abbrev=False,
)

# user show
user_show_parser = user_subparsers.add_parser(
    "show",
    description="Show details of single user.",
    allow_abbrev=False,
)
user_show_parser.add_argument("user_id", help="ID of the user.", type=int)

# user edit
user_edit_parser = user_subparsers.add_parser(
    "edit",
    description=(
        "Edit an existin user and/or their settings "
        "(note: user and settings editing are two independent steps)."
    ),
    allow_abbrev=False,
)
user_edit_parser.add_argument("user_id", help="ID of the user.", type=int)
user_edit_parser.add_argument(
    "--new-email", help="New email address.", required=False
)
user_edit_parser.add_argument(
    "--new-password", help="New password.", required=False
)
user_edit_parser.add_argument(
    "--new-username", help="New user username.", required=False
)
user_edit_parser.add_argument(
    "--new-project-dir",
    help="New value of `project_dir`.",
    required=False,
)
user_edit_parser.add_argument(
    "--new-slurm-user", help="New SLURM username.", required=False
)

user_edit_parser.add_argument(
    "--new-ssh-settings-json",
    help=(
        "Path to JSON file with (a subset of) following settings: "
        "ssh_host, ssh_username, ssh_private_key_path, "
        "ssh_tasks_dir, ssh_jobs_dir."
    ),
    type=str,
    required=False,
)

user_edit_parser_superuser = user_edit_parser.add_mutually_exclusive_group()
user_edit_parser_superuser.add_argument(
    "--make-superuser",
    help="Give superuser privileges to user.",
    action="store_true",
    required=False,
)
user_edit_parser_superuser.add_argument(
    "--remove-superuser",
    help="Remove superuser privileges from user.",
    action="store_true",
    required=False,
)
user_edit_parser_verified = user_edit_parser.add_mutually_exclusive_group()
user_edit_parser_verified.add_argument(
    "--make-verified",
    help="Make user verified.",
    action="store_true",
    required=False,
)
user_edit_parser_verified.add_argument(
    "--remove-verified",
    help="Make user unverified.",
    action="store_true",
    required=False,
)

# user set-groups
user_set_groups_parser = user_subparsers.add_parser(
    "set-groups",
    description=("Reset user-group membership for an existing user."),
    allow_abbrev=False,
)
user_set_groups_parser.add_argument(
    "user_id", help="ID of the user.", type=int
)
user_set_groups_parser.add_argument(
    "group_ids",
    help=(
        "List of the IDs of groups we want the user to be member. "
        "WARNING: this list replaces the current group memberships."
    ),
    type=int,
    nargs="+",
)


# (USER)GROUP GROUP

group_parser = subparsers_main.add_parser(
    "group",
    description="UserGroup commands.",
    allow_abbrev=False,
)
group_subparsers = group_parser.add_subparsers(
    title="Valid sub-commands", dest="subcmd", required=True
)


# group list
group_list_parser = group_subparsers.add_parser(
    "list", description="Get all groups.", allow_abbrev=False
)
group_list_parser.add_argument(
    "--user-ids",
    help="Also return the `user_ids` lists together with the groups",
    action="store_true",
    required=False,
)

# group get
group_get_parser = group_subparsers.add_parser(
    "get", description="Get single group.", allow_abbrev=False
)
group_get_parser.add_argument(
    "group_id", help="ID of the group to get.", type=int
)

# group new
group_new_parser = group_subparsers.add_parser(
    "new", description="Create new group.", allow_abbrev=False
)
group_new_parser.add_argument("name", help="Name of the new group.", type=str)
group_new_parser.add_argument(
    "--viewer-paths",
    help=(
        "List of group's `viewer_paths` (e.g "
        "`--viewer-paths /something /else`)",
    ),
    required=False,
    type=str,
    nargs="+",
)

# group update
group_update_parser = group_subparsers.add_parser(
    "update", description="Update single group.", allow_abbrev=False
)
group_update_parser.add_argument(
    "group_id", help="ID of the group to update.", type=int
)
group_update_parser.add_argument(
    "--new-viewer-paths",
    help=(
        "New list of group `viewer_paths` (e.g "
        "`--new-viewer-paths /something /else`); "
        "note that this replaces the existing one."
    ),
    required=True,
    type=str,
    nargs="+",
    default=None,
)

# group add-user
group_add_user_parser = group_subparsers.add_parser(
    "add-user", description="Add a single user to group.", allow_abbrev=False
)
group_add_user_parser.add_argument(
    "group_id", help="ID of the group to which to add the user.", type=int
)
group_add_user_parser.add_argument(
    "user_id", help="ID of the user to add.", type=int
)

# group remove-user
group_remove_user_parser = group_subparsers.add_parser(
    "remove-user",
    description="Remove a single user from group.",
    allow_abbrev=False,
)
group_remove_user_parser.add_argument(
    "group_id", help="ID of the group to which to remove the user.", type=int
)
group_remove_user_parser.add_argument(
    "user_id", help="ID of the user to remove.", type=int
)
