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
    help="User email address for login (overrides configuration file).",
)
parser_main.add_argument(
    "-p", "--password", help="User password (overrides configuration file)."
)
parser_main.add_argument(
    "--verbose",
    "-v",
    action="store_true",
    default=False,
    help="Change minimal logging level from default to DEBUG.",
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
    "project_id", type=int, help="ID of project to add the new dataset to."
)
project_add_dataset_parser.add_argument(
    "dataset_name", help="Name of new dataset."
)
project_add_dataset_parser.add_argument(
    "--metadata",
    help="Path to file containing dataset metadata in JSON format.",
)
project_add_dataset_parser.add_argument(
    "--type",
    help="Dataset type.",
)
project_add_dataset_parser.add_argument(
    "--make-read-only",
    action="store_true",
    default=False,
    required=False,
    help="Make dataset read-only (not required).",
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
project_edit_parser_read_only = (
    project_edit_parser.add_mutually_exclusive_group()
)
project_edit_parser_read_only.add_argument(
    "--make-read-only",
    help="Set the read-only flag for this project.",
    action="store_true",
    required=False,
)
project_edit_parser_read_only.add_argument(
    "--remove-read-only",
    help="Remove the read-only flag for this project.",
    action="store_true",
    required=False,
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

# dataset add-resource
dataset_add_resource_parser = dataset_subparsers.add_parser(
    "add-resource",
    description="Add resource to existing dataset.",
    allow_abbrev=False,
)
dataset_add_resource_parser.add_argument(
    "project_id", type=int, help="Project ID."
)
dataset_add_resource_parser.add_argument(
    "dataset_id", type=int, help="Dataset ID."
)
dataset_add_resource_parser.add_argument(
    "path", help="Absolute path to resource."
)

# dataset rm-resource
dataset_rm_resource_parser = dataset_subparsers.add_parser(
    "rm-resource",
    description="Remove resource to existing dataset.",
    allow_abbrev=False,
)
dataset_rm_resource_parser.add_argument(
    "project_id", type=int, help="Project ID."
)
dataset_rm_resource_parser.add_argument(
    "dataset_id", type=int, help="Dataset ID."
)
dataset_rm_resource_parser.add_argument(
    "resource_id", type=int, help="Resource ID."
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
dataset_edit_parser.add_argument("--new-type", help="Dataset type.")
dataset_edit_parser.add_argument(
    "--meta-file",
    help="Path to JSON file with new metadata to replace the current ones.",
)
dataset_edit_parser_read = dataset_edit_parser.add_mutually_exclusive_group()
dataset_edit_parser_read.add_argument(
    "--make-read-only",
    help="Set read-only flag of dataset.",
    action="store_true",
    required=False,
)
dataset_edit_parser_read.add_argument(
    "--remove-read-only",
    help="Remove read-only flag of dataset.",
    action="store_true",
    required=False,
)

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

# dataset history
dataset_history_parser = dataset_subparsers.add_parser(
    "history",
    description="Export dataset history as a reproducible worfklow.",
    argument_default=ap.SUPPRESS,
    allow_abbrev=False,
)
dataset_history_parser.add_argument("project_id", type=int, help="Project ID.")
dataset_history_parser.add_argument("dataset_id", type=int, help="Dataset ID.")

# dataset status
dataset_status_parser = dataset_subparsers.add_parser(
    "status",
    description="Extract status of WorkflowTasks associated with dataset.",
    argument_default=ap.SUPPRESS,
    allow_abbrev=False,
)
dataset_status_parser.add_argument("project_id", type=int, help="Project ID.")
dataset_status_parser.add_argument("dataset_id", type=int, help="Dataset ID.")

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

# task check-collection
task_check_collection_parser = task_subparsers.add_parser(
    "check-collection",
    description="Check status of background task collection processes.",
    allow_abbrev=False,
)
task_check_collection_parser.add_argument(
    "state_id",
    help="State ID of the collection (see output of task collect).",
    type=int,
)
task_check_collection_parser.add_argument(
    "--include-logs",
    default=False,
    action="store_true",
    help="Also include task-collection logs.",
)
task_check_collection_parser.add_argument(
    "--do-not-separate-logs",
    dest="do_not_separate_logs",
    help=(
        "Show the task-collection logs in the main output, "
        "instead of a separate field."
    ),
    action="store_true",
    required=False,
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
    "command", help="The command that executes the task."
)
task_new_parser.add_argument("source", help="TBD")
task_new_parser.add_argument(
    "--input-type",
    help="The type of data the task expects as input.",
    default="Any",
)
task_new_parser.add_argument(
    "--output-type",
    help="The type of data the task expects as output.",
    default="Any",
)
task_new_parser.add_argument(
    "--version",
    help="Task version.",
)
task_new_parser.add_argument(
    "--meta-file",
    help="Path to JSON file with additional parameters useful for execution.",
)
task_new_parser.add_argument(
    "--args-schema",
    help="Path to file containing JSON Schema for task arguments.",
)
task_new_parser.add_argument(
    "--args-schema-version",
    help=(
        "Label encoding how the task-arguments JSON Schema was generated "
        "(e.g. `pydantic_v1`)."
    ),
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
task_edit_parser.add_argument("--new-name", help="New task name.")
task_edit_parser.add_argument("--new-command", help="New task command.")
task_edit_parser.add_argument(
    "--new-input-type",
    help="New input type.",
)
task_edit_parser.add_argument(
    "--new-output-type",
    help="New output type.",
)
task_edit_parser.add_argument(
    "--meta-file",
    help=(
        "Path to JSON serialised file containing updates to the current "
        "`meta` dictionary."
    ),
)
task_edit_parser.add_argument(
    "--new-version",
    help="New version.",
)
task_edit_parser.add_argument(
    "--new-args-schema",
    help="Path to file containing the new JSON Schema for task arguments.",
)
task_edit_parser.add_argument(
    "--new-args-schema-version",
    help=(
        "New label encoding how the task-arguments JSON Schema was generated."
    ),
)

# task delete
task_delete_parser = task_subparsers.add_parser(
    "delete",
    description="Delete task.",
    argument_default=ap.SUPPRESS,
    allow_abbrev=False,
)
task_delete_id_or_name_group = task_delete_parser.add_mutually_exclusive_group(
    required=True
)
task_delete_id_or_name_group.add_argument(
    "--id", help="ID of the task to delete.", type=int
)
task_delete_id_or_name_group.add_argument(
    "--name", help="Name of the task to delete."
)
task_delete_parser.add_argument(
    "--version",
    help=(
        "Version of the task to delete "
        "(only accepted in combination with `--name`)."
    ),
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
    "--args-file",
    help=(
        "Path to json serialised file containing the arguments "
        "overrides of the task."
    ),
)
workflow_add_task_parser.add_argument(
    "--meta-file",
    help=(
        "Path to json serialised file containing the meta "
        "overrides of the task."
    ),
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
    "--args-file",
    help=(
        "Path to json serialised file containing the arguments "
        "overrides of the task."
    ),
)
workflow_edit_task_parser.add_argument(
    "--meta-file",
    help=(
        "Path to json serialised file containing the meta "
        "overrides of the task."
    ),
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

# workflow apply
workflow_apply_parser = workflow_subparsers.add_parser(
    "apply",
    description="Apply workflow to dataset.",
    argument_default=ap.SUPPRESS,
    allow_abbrev=False,
)

workflow_apply_parser.add_argument("project_id", type=int)
workflow_apply_parser.add_argument("workflow_id", type=int)
workflow_apply_parser.add_argument("input_dataset_id", type=int)
workflow_apply_parser.add_argument("output_dataset_id", type=int)
workflow_apply_parser.add_argument(
    "--start",
    dest="first_task_index",
    type=int,
    help=(
        "Positional index of the first task to be executed"
        " (starting from 0)."
    ),
    required=False,
)
workflow_apply_parser.add_argument(
    "--end",
    dest="last_task_index",
    type=int,
    help=(
        "Positional index of the last task to be executed"
        " (starting from 0)."
    ),
    required=False,
)
workflow_apply_parser.add_argument(
    "-w",
    "--worker-init",
    help="Command to be run before starting a worker.",
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
    help="Path to a JSON file with the workflow to be imported.",
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
job_show_parser.add_argument(
    "--do-not-separate-logs",
    dest="do_not_separate_logs",
    help="Show the job logs in the main output, instead of a separate field.",
    action="store_true",
    required=False,
)

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

# user register
user_register_parser = user_subparsers.add_parser(
    "register",
    description="Register a new user with the Fractal server.",
    allow_abbrev=False,
)
user_register_parser.add_argument(
    "new_email", help="Email to be used as username."
)
user_register_parser.add_argument(
    "new_password", help="Password for the new user."
)
user_register_parser.add_argument(
    "--cache-dir",
    help=(
        "User's cache directory absolute path "
        "(necessary for workflow execution when using the SLURM backend)."
    ),
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
    description="Edit details of single user.",
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
    "--new-cache-dir",
    help=(
        "New user's cache directory absolute path "
        "(necessary for workflow execution when using the SLURM backend)."
    ),
    required=False,
)
user_edit_parser.add_argument(
    "--new-slurm-user", help="New SLURM username.", required=False
)
user_edit_parser.add_argument(
    "--new-username", help="New user username.", required=False
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

# user delete
user_delete_parser = user_subparsers.add_parser(
    "delete",
    description="Delete a single user.",
    allow_abbrev=False,
)
user_delete_parser.add_argument("user_id", help="ID of the user.", type=int)
