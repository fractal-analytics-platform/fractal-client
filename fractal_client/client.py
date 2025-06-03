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

from httpx import ConnectError

from . import cmd
from .authclient import AuthClient
from .authclient import AuthenticationError
from .config import settings
from .interface import Interface
from .parser import parser_main


def _verify_authentication_branch(
    *,
    username: str | None,
    password: str | None,
    token_path: str | None,
) -> None:
    """
    Fail if credentials are not either username&password or token

    Arguments:
        username: Username
        password: Password
        token_path: Path of token
    """
    which_parameters_are_set = (
        bool(username),
        bool(password),
        bool(token_path),
    )
    valid_cases = (
        (True, True, False),
        (False, False, True),
    )
    if which_parameters_are_set not in valid_cases:
        msg = (
            "Invalid authentication credentials. "
            "You should either set username&password or the token path.\n\n"
            "You can set these variables in multiple ways "
            "(see `fractal --help`):\n"
            "  1. Through command-line arguments.\n"
            "  2. Through environment variables.\n"
            "  3. Through environment variables in a `.fractal.env` file.\n"
        )
        raise ValueError(msg)


def handle(cli_args: list[str] = sys.argv) -> Interface:

    args = parser_main.parse_args(cli_args[1:])

    # Set logging level
    if args.debug is True:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    show_args = "\n".join(
        [
            f"    {k}: {v}"
            if not (k == "password" and v is not None)
            else "    password: ***"
            for k, v in args.__dict__.items()
        ]
    )
    logging.debug(f"\nArguments:\n{show_args}")

    if args.cmd is not None:
        handler = getattr(cmd, args.cmd)
    else:
        parser_main.print_help()
        sys.exit(1)

    try:
        # Make a copy of vars(args), and remove cmd (which is not a relevant
        # argument for functions called with **kwargs)
        kwargs = vars(args).copy()
        kwargs.pop("cmd")
        fractal_server = (
            kwargs.pop("fractal_server") or settings.FRACTAL_SERVER
        )
        logging.debug(f"Fractal server URL: {fractal_server}")
        if fractal_server is None:
            return Interface(
                data=(
                    "Missing argument: You should set the "
                    "fractal-server URL (see `fractal --help`)."
                ),
                retcode=1,
            )

        if args.cmd == "version":
            interface = handler(fractal_server, **kwargs)
        else:
            # Extract (and remove) credentials-related variables from kwargs
            username = kwargs.pop("user") or settings.FRACTAL_USER
            password = kwargs.pop("password") or settings.FRACTAL_PASSWORD
            token_path = (
                kwargs.pop("token_path") or settings.FRACTAL_TOKEN_PATH
            )
            try:
                _verify_authentication_branch(
                    username=username,
                    password=password,
                    token_path=token_path,
                )
            except ValueError as e:
                return Interface(data=str(e), retcode=1)
            # Read token from file
            if token_path is not None:
                try:
                    with open(token_path) as f:
                        token = f.read().strip("\n")
                except Exception as e:
                    msg = (
                        f"Could not read token from {token_path=}.\n"
                        f"Original error:\n{str(e)}"
                    )
                    return Interface(data=msg, retcode=1)

            else:
                token = None

            with AuthClient(
                fractal_server=fractal_server,
                username=username,
                password=password,
                token=token,
            ) as client:
                interface = handler(client, **kwargs)
    except AuthenticationError as e:
        return Interface(retcode=1, data=e.args[0])
    except ConnectError as e:
        return Interface(
            retcode=1,
            data=(
                f"ConnectError at {e.request.url}\n"
                f"Original error: '{e.args[0]}'\n"
                f"Hint: is {settings.FRACTAL_SERVER} alive?"
            ),
        )
    except Exception as e:
        return Interface(retcode=1, data=str(e))

    return interface


def main():
    interface = handle()
    interface.show()
    sys.exit(interface.retcode)


if __name__ == "__main__":
    main()
