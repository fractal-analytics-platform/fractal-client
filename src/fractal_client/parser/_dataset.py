import argparse as ap


def add_dataset_parser(subparsers):
    dataset_parser = subparsers.add_parser(
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
    dataset_edit_parser.add_argument(
        "--new-name", type=str, required=True, help="New name of dataset."
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
