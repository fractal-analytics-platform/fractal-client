import argparse as ap


def add_workflow_parser(subparsers):
    workflow_parser = subparsers.add_parser(
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
    workflow_delete_parser.add_argument("workflow_id", type=int, help="Workflow ID.")

    # workflow add task
    workflow_add_task_parser = workflow_subparsers.add_parser(
        "add-task",
        description="Append a single task to the task list of a workflow.",
        allow_abbrev=False,
    )
    workflow_add_task_parser.add_argument("project_id", type=int, help="Project ID.")
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
            "Version of task to add (only accepted in combination with --task-name)."
        ),
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
    workflow_edit_task_parser.add_argument("project_id", type=int, help="Project ID.")
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
        help=("Path to JSON file containing the type filters of this workflow task."),
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
    workflow_remove_task_parser.add_argument("project_id", type=int, help="Project ID.")
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

    # workflow import-from-template
    workflow_import_from_template_parser = workflow_subparsers.add_parser(
        "import-from-template",
        description="Import workflow to project from template.",
        allow_abbrev=False,
    )
    workflow_import_from_template_parser.add_argument(
        "project_id",
        type=int,
        help="ID of the project where the workflow will be imported.",
    )
    workflow_import_from_template_parser.add_argument(
        "template_id",
        type=int,
        help="ID of the template from which the workflow will be imported.",
    )
    workflow_import_from_template_parser.add_argument(
        "--name",
        type=str,
        help=("Name of the new workflow (if set, overrides the one in the template)."),
        required=False,
    )
