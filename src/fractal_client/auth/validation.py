import argparse as ap
import sys
from collections.abc import Callable

import fractal_client.cmd as cmd_module
from fractal_client.config import settings

from .info import AuthInfo

ERROR_MSG = (
    "Invalid authentication credentials. "
    "You should either set username&password or the token path.\n\n"
    "You can set these variables in multiple ways "
    "(see `fractal --help`):\n"
    "  1. Through command-line arguments (e.g. `-u`, `-p`, `--token-path`).\n"
    "  2. Through environment variables (e.g. `FRACTAL_USER`).\n"
    "  3. Through environment variables in a `.fractal.env` file.\n"
)


def get_fractal_server(args: ap.Namespace) -> str:
    fractal_server: str | None = args.fractal_server or settings.FRACTAL_SERVER
    if fractal_server is None:
        print(
            "Missing argument: You should set the "
            "fractal-server URL (see `fractal --help`)."
        )
        sys.exit(1)
    else:
        fractal_server = fractal_server.rstrip("/")
        return fractal_server


def get_cmd_handler(*, parser: ap.ArgumentParser, args: ap.Namespace) -> Callable:
    cmd: str | None = args.cmd
    if cmd is None:
        parser.print_help()
        sys.exit(1)
    elif (cmd_handler := getattr(cmd_module, cmd, None)) is not None:
        return cmd_handler
    else:
        print(f"An internal error took place while processing {cmd=}.")
        sys.exit(1)


def get_auth_info(*, parser: ap.ArgumentParser, args: ap.Namespace) -> AuthInfo:
    user: str | None = args.user or settings.FRACTAL_SERVER
    password: str | None = args.password or settings.FRACTAL_PASSWORD
    token_path: str | None = args.token_path or settings.FRACTAL_TOKEN_PATH

    if bool(user) != bool(password):
        print(ERROR_MSG)
        sys.exit(1)

    if bool(user) and bool(token_path):
        print(ERROR_MSG)
        sys.exit(1)

    return AuthInfo(
        user=user,
        password=password,
        token_path=token_path,
    )
