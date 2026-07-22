def add_template_parser(subparsers):
    template_parser = subparsers.add_parser(
        "template",
        description="Template commands.",
        allow_abbrev=False,
    )
    template_subparsers = template_parser.add_subparsers(
        title="Valid sub-commands", dest="subcmd", required=True
    )

    # template show
    template_show_parser = template_subparsers.add_parser(
        "show",
        description="Show single template.",
        allow_abbrev=False,
    )
    template_show_parser.add_argument(
        "template_id", help="ID of the template to show.", type=int
    )

    # template new
    template_new_parser = template_subparsers.add_parser(
        "new",
        description="Create new template.",
        allow_abbrev=False,
    )

    template_new_from_workflow_or_import = (
        template_new_parser.add_mutually_exclusive_group(required=True)
    )
    template_new_from_workflow_or_import.add_argument(
        "--workflow-id",
        help="ID of the workflow from which the new template will be built.",
        type=int,
    )
    template_new_from_workflow_or_import.add_argument(
        "--json-file",
        help="Path to a JSON file with the template to be imported.",
    )

    template_new_parser.add_argument(
        "--name", help="New template name.", required=False
    )
    template_new_parser.add_argument(
        "--version", help="New template version.", required=False
    )
    template_new_parser.add_argument(
        "--user-group-id",
        help=(
            "ID of user group which should be granted access to the new template."
        ),
        required=False,
    )

    # template delete
    template_delete_parser = template_subparsers.add_parser(
        "delete",
        description="Delete single template.",
        allow_abbrev=False,
    )
    template_delete_parser.add_argument(
        "template_id", help="ID of the template to delete.", type=int
    )

    # template export
    template_export_parser = template_subparsers.add_parser(
        "export",
        description="Export single template.",
        allow_abbrev=False,
    )
    template_export_parser.add_argument(
        "template_id",
        help="ID of the template to export.",
    )
    template_export_parser.add_argument(
        "json_file",
        help="Path where to export the template.",
    )
