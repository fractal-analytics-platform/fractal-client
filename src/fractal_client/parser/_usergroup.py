def add_usergroup_parser(subparsers):
    group_parser = subparsers.add_parser(
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
    group_new_parser.add_argument(
        "name", help="Name of the new group.", type=str
    )

    # group add-user
    group_add_user_parser = group_subparsers.add_parser(
        "add-user",
        description="Add a single user to group.",
        allow_abbrev=False,
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
        "group_id",
        help="ID of the group to which to remove the user.",
        type=int,
    )
    group_remove_user_parser.add_argument(
        "user_id", help="ID of the user to remove.", type=int
    )
