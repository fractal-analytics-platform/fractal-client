"""
Copyright 2022 (C) Friedrich Miescher Institute for Biomedical Research and
University of Zurich

Original authors:
Jacopo Nespolo <jacopo.nespolo@exact-lab.it>
Marco Franzon <marco.franzon@exact-lab.it>
Tommaso Comparin <tommaso.comparin@exact-lab.it>

This file is part of Fractal and was originally developed by eXact lab S.r.l.
<exact-lab.it> under contract with Liberali Lab from the Friedrich Miescher
Institute for Biomedical Research and Pelkmans Lab from the University of
Zurich.
"""
import argparse as ap


parser_main = ap.ArgumentParser(description="Fractal Analytics Framework")

parser_main.add_argument(
    "-u",
    "--user",
    help="User email address for login (overrides configuration file)",
)
parser_main.add_argument(
    "-p", "--password", help="User password (overrides configuration file)"
)
parser_main.add_argument("-s", "--slurm_user", help="Slurm user")
parser_main.add_argument(
    "-v",
    action="count",
    default=0,
    help="Use one or more time to increase verbosity level",
)
parser_main.add_argument(
    "--batch",
    default=False,
    action="store_true",
    help=(
        "Return output suitable for scripting, e.g., "
        "only the id of items created instead of the full object."
    ),
)
parser_main.add_argument(
    "-j",
    "--json",
    default=False,
    action="store_true",
    help="Output raw json",
)

subparsers_main = parser_main.add_subparsers(title="Commands:", dest="cmd")

# REGISTER GROUP
register_parser = subparsers_main.add_parser(
    "register", help="Register with the Fractal server"
)
register_parser.add_argument("email", help="Email to be used as username")
register_parser.add_argument(
    "slurm_user", help="Username to login into Slurm cluster"
)
register_parser.add_argument(
    "-p",
    "--password",
    help=("Password for the new user"),
)
# PROJECT GROUP
project_parser = subparsers_main.add_parser("project", help="project commands")
project_subparsers = project_parser.add_subparsers(
    title="Valid subcommands:", dest="subcmd", required=True
)

# project new
project_new_parser = project_subparsers.add_parser(
    "new", help="Create new project"
)
project_new_parser.add_argument("name", help="Name of new project", type=str)
project_new_parser.add_argument(
    "path",
    help=(
        "Project directory of new project. "
        "New datasets will be written here."
    ),
)
project_new_parser.add_argument(
    "-d",
    "--dataset",
    help=(
        "Name of new dataset to create. "
        "By default, the dataset `default` is created."
    ),
)

# project list
project_list_parser = project_subparsers.add_parser(
    "list", help="List projects"
)

# project show
project_show_parser = project_subparsers.add_parser(
    "show", help="Show details of a single project"
)
project_show_parser.add_argument(
    "project_id", type=int, help="Id of project to show"
)


# project delete
project_delete_parser = project_subparsers.add_parser(
    "delete", help="Delete project"
)

# project add-dataset
project_add_dataset_parser = project_subparsers.add_parser(
    "add-dataset", help="Add dataset to project"
)
project_add_dataset_parser.add_argument(
    "project_id", type=int, help="Id of project to add the new dataset to"
)
project_add_dataset_parser.add_argument(
    "dataset_name", help="Name of new dataset"
)
project_add_dataset_parser.add_argument(
    "--metadata",
    help="Path to file containing dataset metadata in JSON format.",
)


# DATASET GROUP
dataset_parser = subparsers_main.add_parser("dataset", help="dataset commands")
dataset_subparsers = dataset_parser.add_subparsers(
    title="Valid subcommands:", dest="subcmd", required=True
)

# dataset add-resource
dataset_add_resource_parser = dataset_subparsers.add_parser(
    "add-resource", help="Add resource to existing dataset"
)
dataset_add_resource_parser.add_argument(
    "project_id", type=int, help="Project id"
)
dataset_add_resource_parser.add_argument(
    "dataset_id", type=int, help="Dataset id"
)
dataset_add_resource_parser.add_argument("path", help="Path to resource")
dataset_add_resource_parser.add_argument(
    "-g", "--glob-pattern", help="Glob pattern"
)

# dataset add-resource
dataset_rm_resource_parser = dataset_subparsers.add_parser(
    "rm-resource", help="Remove resource to existing dataset"
)
dataset_rm_resource_parser.add_argument(
    "project_id", type=int, help="Project id"
)
dataset_rm_resource_parser.add_argument(
    "dataset_id", type=int, help="Dataset id"
)
dataset_rm_resource_parser.add_argument(
    "resource_id", type=int, help="Resource id"
)

# dataset edit
dataset_edit_parser = dataset_subparsers.add_parser(
    "edit", help="Edit dataset", argument_default=ap.SUPPRESS
)
dataset_edit_parser.add_argument("project_id", type=int, help="Project id")
dataset_edit_parser.add_argument("dataset_id", type=int, help="Dataset id")
dataset_edit_parser.add_argument("--name", help="New name of dataset")
dataset_edit_parser.add_argument("--path", help="New path of dataset")
dataset_edit_parser.add_argument(
    "--metadata",
    help=(
        "Path to file containing dataset metadata in JSON format. "
        "(Set to `none` to clear)"
    ),
)
dataset_edit_parser.add_argument(
    "--read-only",
    dest="read_only",
    help="Set read-only flag of dataset",
    action="store_true",
)
dataset_edit_parser.add_argument(
    "--read-write",
    dest="read_only",
    help="Set read-only flag of dataset (0 for False, 1 for True)",
    action="store_false",
)
dataset_edit_parser.add_argument("-t", "--type", help="Dataset type")

# dataset show
dataset_show_parser = dataset_subparsers.add_parser(
    "show", help="Show dataset", argument_default=ap.SUPPRESS
)
dataset_show_parser.add_argument("project_id", type=int, help="Project id")
dataset_show_parser.add_argument("dataset_id", type=int, help="Dataset id")


# TASK GROUP
task_parser = subparsers_main.add_parser("task", help="task commands")
task_subparsers = task_parser.add_subparsers(
    title="Valid subcommands:", dest="subcmd", required=True
)

# task list
task_list_parser = task_subparsers.add_parser("list", help="List tasks")

# task collect
task_collect_parser = task_subparsers.add_parser(
    "collect",
    help="Install and collect all tasks a pip installable package exposes",
)
task_collect_parser.add_argument(
    "package",
    help="Package name or path to local package",
)
task_collect_parser.add_argument(
    "--python-version",
    help="Select the python version to use for this package",
)
task_collect_parser.add_argument(
    "--package-version",
    help="Select the package version",
)
task_collect_parser.add_argument(
    "--package-extras",
    help=(
        "Comma separated list of extra components for the package to be "
        "installed, e.g., `collect fractal-tasks-core "
        "--package-extras=torch,tensorflow` will trigger the installation of "
        "`fractal-tasks-core[torch,tensorflow]`"
    ),
)
task_collect_parser.add_argument(
    "--private",
    default=False,
    action="store_true",
    help="Intall tasks as private to the user (as opposed to global)",
)

# task check-collection
task_check_collection_parser = task_subparsers.add_parser(
    "check-collection",
    help="Check status of background task collection processes",
)
task_check_collection_parser.add_argument(
    "state_id",
    help="State ID of the collection (see output of task collect)",
)
task_check_collection_parser.add_argument(
    "--verbose",
    default=False,
    action="store_true",
    help="Output more verbose output",
)


# task edit
task_edit_parser = task_subparsers.add_parser(
    "edit", help="Edit task", argument_default=ap.SUPPRESS
)
task_edit_parser.add_argument(
    "task_id_or_name", help="ID or name of task to edit", type=str
)
task_edit_parser.add_argument("--name", help="New task name")
task_edit_parser.add_argument(
    "--resource-type",
    choices=["task", "workflow"],
    help="New resource type",
)
task_edit_parser.add_argument(
    "--input-type",
    help="New input type",
)
task_edit_parser.add_argument(
    "--output-type",
    help="New resource type",
)
task_edit_parser.add_argument(
    "--default-args", help="Filename containing JSON encoded default arguments"
)

# WORKFLOW GROUP

# workflow new
workflow_parser = subparsers_main.add_parser(
    "workflow", help="workflow commands"
)
workflow_subparsers = workflow_parser.add_subparsers(
    title="Valid subcommand", dest="subcmd", required=True
)
workflow_new_parser = workflow_subparsers.add_parser(
    "new", help="Create new workflow"
)
workflow_new_parser.add_argument(
    "name",
    help="Workflow name (must be unique, and not only made of numbers only)",
)
workflow_new_parser.add_argument(
    "project_id",
    help="Project ID",
)

# workflow list
workflow_list_parser = workflow_subparsers.add_parser(
    "list", help="List workflows for given project"
)
workflow_list_parser.add_argument(
    "project_id",
    help="Project ID",
)

# workflow delete
workflow_new_parser = workflow_subparsers.add_parser(
    "delete", help="Delete workflow"
)
workflow_new_parser.add_argument(
    "id",
    help="Workflow id",
)

# workflow show
workflow_new_parser = workflow_subparsers.add_parser(
    "show", help="Show a workflow"
)
workflow_new_parser.add_argument(
    "id",
    help="Workflow id",
)

# workflow add task
workflow_add_task_parser = workflow_subparsers.add_parser(
    "add-task", help="Add a new task to a specific workflow"
)
workflow_add_task_parser.add_argument(
    "id",
    help="Workflow id",
)
workflow_add_task_parser.add_argument(
    "task_id_or_name", help="ID or name of the new task", type=str
)
workflow_add_task_parser.add_argument(
    "--order", help="Order of this task within the workflow's task list"
)
workflow_add_task_parser.add_argument(
    "--args-file",
    help=(
        "Path to a json serialised file containing the arguments "
        "ovverrides of the task"
    ),
)
workflow_add_task_parser.add_argument(
    "--meta-file",
    help=(
        "Path to a json serialised file containing the meta"
        "ovverrides of the task"
    ),
)


# workflow edit task
workflow_edit_task_parser = workflow_subparsers.add_parser(
    "edit-task", help="Edit a task within a specific workflow"
)
workflow_edit_task_parser.add_argument(
    "id",
    help="Workflow id",
)
workflow_edit_task_parser.add_argument(
    "workflow_task_id",
    help="Workflow task Id, the Id of a task inside the list of tasks",
)
workflow_edit_task_parser.add_argument(
    "--args-file",
    help=(
        "Path to a json serialised file containing the arguments "
        "ovverrides of the task"
    ),
)
workflow_edit_task_parser.add_argument(
    "--meta-file",
    help=(
        "Path to a json serialised file containing the meta"
        "ovverrides of the task"
    ),
)


# workflow remove task
workflow_remove_task_parser = workflow_subparsers.add_parser(
    "rm-task", help="Remove a task in a specific workflow"
)
workflow_remove_task_parser.add_argument(
    "id",
    help="Workflow id",
)
workflow_remove_task_parser.add_argument(
    "workflow_task_id",
    help="Workflow task Id, the Id of a task inside the list of tasks",
)

# workflow edit
workflow_edit_parser = workflow_subparsers.add_parser(
    "edit", help="Edit workflow", argument_default=ap.SUPPRESS
)
workflow_edit_parser.add_argument(
    "id",
    help="Workflow id",
)
workflow_edit_parser.add_argument("--name", help="New workflow name")

workflow_edit_parser.add_argument(
    "--project-id",
    help="Change the project associated with the current workflow",
)

# workflow apply
workflow_apply_parser = workflow_subparsers.add_parser(
    "apply", help="Apply workflow to dataset", argument_default=ap.SUPPRESS
)
workflow_apply_parser.add_argument("workflow_id")
workflow_apply_parser.add_argument("input_dataset_id")
workflow_apply_parser.add_argument(
    "-o", "--output_dataset_id", help="Output dataset id"
)
workflow_apply_parser.add_argument(
    "--overwrite-input",
    default=False,
    action="store_true",
    help="Allow overwriting the content of the input dataset",
)
workflow_apply_parser.add_argument(
    "-p",
    "--project-id",
    help="Id of project the workflow and dataset belong to",
)
workflow_apply_parser.add_argument(
    "-w",
    "--worker-init",
    help="Command to be run before starting a worker",
)


# JOB GROUP

job_parser = subparsers_main.add_parser("job", help="job commands")
job_subparsers = job_parser.add_subparsers(
    title="Valid subcommand", dest="subcmd", required=True
)

# job list
job_list_parser = job_subparsers.add_parser(
    "list", help="List jobs for given project"
)
job_list_parser.add_argument(
    "project_id",
    help="Project ID",
)

# job status
job_status_parser = job_subparsers.add_parser(
    "status",
    help="Query status of workflow-execution job",
    argument_default=ap.SUPPRESS,
)
job_status_parser.add_argument(
    "job_id",
    help="Id of the job",
)
job_status_parser.add_argument(
    "--do-not-separate-logs",
    dest="do_not_separate_logs",
    help="Show the job logs in the main output, instead of a separate field",
    action="store_true",
)

# job download-logs
job_download_logs_parser = job_subparsers.add_parser(
    "download-logs",
    help="Download full folder of workflow-execution job",
)
job_download_logs_parser.add_argument(
    "job_id",
    help="Id of the job",
)
job_download_logs_parser.add_argument(
    "--output",
    help="Path of the output folder",
)


# VERSION GROUP
version_parser = subparsers_main.add_parser(
    "version", help="Print verison and exit"
)
