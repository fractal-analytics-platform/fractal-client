import argparse as ap
import sys
from collections.abc import Callable

import fractal_client.cmd as cmd_module
from fractal_client.config import settings


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
