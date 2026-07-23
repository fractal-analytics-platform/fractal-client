def add_user_parser(subparsers):
    user_parser = subparsers.add_parser(
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
    user_register_parser.add_argument("new_email", help="Email to be used as username.")
    user_register_parser.add_argument("new_password", help="Password for the new user.")
    user_register_parser.add_argument(
        "new_project_dir",
        help="User-writeable base folder, used e.g. for default `zarr_dir` paths.",
    )

    user_register_parser.add_argument(
        "--superuser",
        help="Give superuser privileges to the new user.",
        action="store_true",
        required=False,
    )

    # user list
    user_subparsers.add_parser(
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
            "Edit an existing user and/or their settings "
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
        "--add-project-dir",
        help="New folder to add to user `project_dirs`.",
        required=False,
    )
    user_edit_parser.add_argument(
        "--remove-project-dir",
        help="Folder to remove from user `project_dirs`.",
        required=False,
    )
    user_edit_parser.add_argument(
        "--new-profile-id", help="New value of `profile_id`", required=False
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
    user_set_groups_parser.add_argument("user_id", help="ID of the user.", type=int)
    user_set_groups_parser.add_argument(
        "group_ids",
        help=(
            "List of the IDs of groups we want the user to be member. "
            "WARNING: this list replaces the current group memberships."
        ),
        type=int,
        nargs="+",
    )
