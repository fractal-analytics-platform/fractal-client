def add_auth_parser(subparsers):
    auth_parser = subparsers.add_parser(
        "auth",
        description="Authentication commands.",
        allow_abbrev=False,
    )
    auth_subparsers = auth_parser.add_subparsers(
        title="Valid sub-commands", dest="subcmd", required=True
    )

    check_token_parser = auth_subparsers.add_parser(  # noqa: F841
        "check-token",
        description="Check whether a valid token is available.",
        allow_abbrev=False,
    )

    set_token_parser = auth_subparsers.add_parser(  # noqa: F841
        "set-token",
        description=(
            "Write a valid token to disk (using the default path or "
            "the user-provided `--token-path` one)."
        ),
        allow_abbrev=False,
    )

    clear_token_parser = auth_subparsers.add_parser(  # noqa: F841
        "clear-token",
        description="Remove the file storing the token, if valid.",
        allow_abbrev=False,
    )
