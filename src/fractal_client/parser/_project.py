def add_project_parser(subparsers):
    # PROJECT GROUP
    project_parser = subparsers.add_parser(
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
    project_subparsers.add_parser(
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
        "--project-dir",
        type=str,
        help=(
            "Choose which project directory your dataset zarr directory is placed "
            "in. To add additional project directory choices, contact an admin."
        ),
        required=False,
    )
    project_add_dataset_parser.add_argument(
        "--zarr-subfolder",
        type=str,
        help=(
            "Specify where in your project directory the dataset zarr directory "
            "should be. This is a path relative to the project directory and "
            "needs to stay within the chosen project directory. By default, "
            "Fractal will create a folder for the project with a subfolder for "
            "the dataset."
        ),
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
