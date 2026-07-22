import argparse as ap


def get_main_parser() -> ap.ArgumentParser:
    parser_main = ap.ArgumentParser(
        description="Command-line interface for Fractal Client.",
        allow_abbrev=False,
    )

    parser_main.add_argument(
        "-u",
        "--user",
        help=(
            "User email address for login (overrides `FRACTAL_USER` "
            "environment variable)."
        ),
    )
    parser_main.add_argument(
        "-p",
        "--password",
        help="User password (overrides `FRACTAL_PASSWORD` environment variable).",
    )

    parser_main.add_argument(
        "--token-path",
        help="User token (overrides `FRACTAL_TOKEN_PATH` environment variable).",
    )
    parser_main.add_argument(
        "--fractal-server",
        help=(
            "URL of Fractal server (overrides `FRACTAL_SERVER` environment variable). "
            "Example: https://fractal.example.org/backend"
        ),
    )

    parser_main.add_argument(
        "--fractal-web",
        help=(
            "URL of Fractal web (overrides `FRACTAL_WEB` environment variable)."
            "Example: https://fractal.example.org"
        ),
    )

    parser_main.add_argument(
        "--debug",
        action="store_true",
        default=False,
        help="Change minimal logging level from INFO to DEBUG.",
    )
    parser_main.add_argument(
        "--batch",
        default=False,
        action="store_true",
        help=(
            "Return output suitable for scripting, e.g., "
            "only the ID of items created instead of the full object."
        ),
    )

    return parser_main
