def add_version_parser(subparsers):
    subparsers.add_parser(
        "version",
        description="Print version and exit.",
        allow_abbrev=False,
    )
