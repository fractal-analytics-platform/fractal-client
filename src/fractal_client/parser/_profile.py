def add_profile_parser(subparsers):
    profile_parser = subparsers.add_parser(
        "profile",
        description="Profile commands.",
        allow_abbrev=False,
    )
    profile_subparsers = profile_parser.add_subparsers(
        title="Valid sub-commands", dest="subcmd", required=True
    )

    # profile create
    profile_new_parser = profile_subparsers.add_parser(
        "new", description="Create new profile.", allow_abbrev=False
    )
    profile_new_parser.add_argument("resource_id", help="TBD")
    profile_new_parser.add_argument("json_file", help="TBD")
