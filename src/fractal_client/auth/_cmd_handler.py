"""
Get a Python function based on a CLI subcommand.
"""

import argparse as ap
import sys
from collections.abc import Callable

import fractal_client.cmd as cmd_module


def get_cmd_handler(*, parser_help: str, args: ap.Namespace) -> Callable:
    cmd: str | None = args.cmd
    if cmd is None:
        sys.exit(parser_help)
    elif (cmd_handler := getattr(cmd_module, cmd, None)) is not None:
        return cmd_handler
    else:
        sys.exit(f"An internal error took place while processing {cmd=}.")
