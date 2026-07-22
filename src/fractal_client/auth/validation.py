import argparse as ap
import sys
from collections.abc import Callable

import fractal_client.cmd as cmd_module
from fractal_client.config import settings


def get_fractal_server(args: ap.Namespace) -> str:
    fractal_server: str | None = args.fractal_server or settings.FRACTAL_SERVER
    if fractal_server is None:
        sys.exit(
            "Missing argument: You should set the "
            "fractal-server URL (see `fractal --help`)."
        )
    else:
        fractal_server = fractal_server.rstrip("/")
        return fractal_server


def get_cmd_handler(*, parser: ap.ArgumentParser, args: ap.Namespace) -> Callable:
    cmd: str | None = args.cmd
    if cmd is None:
        sys.exit(parser.format_help())
    elif (cmd_handler := getattr(cmd_module, cmd, None)) is not None:
        return cmd_handler
    else:
        sys.exit(f"An internal error took place while processing {cmd=}.")
