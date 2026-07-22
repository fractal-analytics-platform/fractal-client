def add_resource_parser(subparsers):
    resource_parser = subparsers.add_parser(
        "resource",
        description="Resource commands.",
        allow_abbrev=False,
    )
    resource_subparsers = resource_parser.add_subparsers(
        title="Valid sub-commands", dest="subcmd", required=True
    )

    # resource create
    resource_new_parser = resource_subparsers.add_parser(
        "new", description="Create new resource.", allow_abbrev=False
    )
    resource_new_parser.add_argument("json_file", help="TBD")
