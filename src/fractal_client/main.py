"""
Copyright 2022 (C) Friedrich Miescher Institute for Biomedical Research and
University of Zurich

Original authors:
Jacopo Nespolo <jacopo.nespolo@exact-lab.it>
Marco Franzon <marco.franzon@exact-lab.it>
Tommaso Comparin <tommaso.comparin@exact-lab.it>

This file is part of Fractal and was originally developed by eXact lab S.r.l.
<exact-lab.it> under contract with Liberali Lab from the Friedrich Miescher
Institute for Biomedical Research and Pelkmans Lab from the University of
Zurich.
"""

import logging
import sys
from collections.abc import Callable

from httpx2 import ConnectError

from .auth._args import AuthInfo
from .auth._args import get_auth_info
from .auth._args import get_fractal_server
from .auth._client import AuthClient
from .auth._client import AuthenticationError
from .auth._cmd_handler import get_cmd_handler
from .config import settings
from .interface import Interface
from .parser import parser_main


def handle(cli_args: list[str]) -> Interface:
    # Parse CLI arguments
    args = parser_main.parse_args(cli_args[1:])

    # Set logging level
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    # Validate backend URL
    fractal_server = get_fractal_server(args)
    logging.debug(f"Fractal server URL: {fractal_server}")

    # Print debug info
    show_args = "\n".join(
        [
            f"    {k}: {v}"
            if not (k == "password" and v is not None)
            else "    password: ****"
            for k, v in args.__dict__.items()
        ]
    )
    logging.debug(f"\nArguments:\n{show_args}")

    # Get command handler function
    cmd_handler: Callable = get_cmd_handler(
        args=args,
        parser_help=parser_main.format_help(),
    )

    try:
        if args.cmd == "version":
            interface = cmd_handler(fractal_server)
        else:
            auth_info: AuthInfo = get_auth_info(args=args)
            with AuthClient(
                fractal_server=fractal_server, auth_info=auth_info
            ) as auth_client:
                # Make a copy of vars(args), and remove `cmd` (which is not a relevant
                # argument for functions called with **kwargs)
                kwargs = vars(args).copy()
                for key in ("cmd", "user", "password", "token_path", "fractal_server"):
                    kwargs.pop(key)
                interface = cmd_handler(auth_client, **kwargs)
    except ConnectError as e:
        return Interface(
            retcode=1,
            data=(
                f"ConnectError at {e.request.url}\n"
                f"Original error: '{e.args[0]}'\n"
                f"Hint: is {settings.FRACTAL_SERVER} alive?"
            ),
        )
    except AuthenticationError as e:
        return Interface(retcode=1, data=e.args[0])
    except Exception as e:
        return Interface(retcode=1, data=str(e))

    return interface


def main():
    interface = handle(sys.argv)
    interface.show()
    sys.exit(interface.retcode)
