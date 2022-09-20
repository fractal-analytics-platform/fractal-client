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
    "-p", "--password", help="User password (overrides cnofiguration file)"
)
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

# PROJECT GROUP
project_parser = subparsers_main.add_parser("project", help="project commands")
project_subparsers = project_parser.add_subparsers(
    title="Valid subcommands:", dest="subcmd", required=True
)

# project new
project_new_parser = project_subparsers.add_parser(
    "new", help="Create new project"
)
project_new_parser.add_argument("name", help="Name of new project")
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

# task new
task_new_parser = task_subparsers.add_parser("new", help="List tasks")
task_new_parser.add_argument("name", help="Task name (must be unique)")
task_new_parser.add_argument(
    "resource_type", choices=["task", "workflow"], help="Resource type"
)
task_new_parser.add_argument("input_type", help="Dataset input type")
task_new_parser.add_argument("output_type", help="Dataset output type")
task_new_parser.add_argument(
    "--module",
    help="Module path, e.g., `module.submodule:task_function`",
)
task_new_parser.add_argument(
    "--default_args", help="Default arguments, as a JSON-encoded string"
)
task_new_parser.add_argument("--subtask-list", help="Subtask list")


# task edit
task_edit_parser = task_subparsers.add_parser(
    "edit", help="Edit task", argument_default=ap.SUPPRESS
)
task_edit_parser.add_argument("task_id", help="ID of task to edit")
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

# task add-subtask
task_add_subtask_parser = task_subparsers.add_parser(
    "add-subtask", help="Edit task", argument_default=ap.SUPPRESS
)
task_add_subtask_parser.add_argument(
    "parent_task_id", help="ID of task to which the subtask will be added"
)
task_add_subtask_parser.add_argument(
    "subtask_id", help="ID of task to add as a subtask"
)
task_add_subtask_parser.add_argument(
    "--args-file",
    help="Path to file containing JSON serialised task arguments",
)


# task apply
task_apply_parser = task_subparsers.add_parser(
    "apply", help="Apply task to a dataset"
)
task_apply_parser.add_argument(
    "project_id", help="ID of project on which to apply task"
)
task_apply_parser.add_argument("input_dataset_id", help="ID of input dataset")
task_apply_parser.add_argument(
    "output_dataset_id", help="ID of output dataset"
)
task_apply_parser.add_argument(
    "workflow_id", help="ID of taks/workflow to apply"
)
task_apply_parser.add_argument(
    "--overwrite_input",
    action="store_true",
    default=False,
    help="Overwrite the input dataset",
)

# VERSION GROUP
version_parser = subparsers_main.add_parser(
    "version", help="Print verison and exit"
)


if __name__ == "__main__":
    args = parser_main.parse_args()

    from devtools import debug

    debug(args)
