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
from sys import argv

from httpx import ConnectError

from . import cmd
from .authclient import AuthClient
from .authclient import AuthenticationError
from .config import settings
from .interface import Interface
from .parser import parser_main


class MissingCredentialsError(RuntimeError):
    pass


def _validate_credentials(
    *,
    username: str | None,
    password: str | None,
    token_path: str | None,
) -> None:
    """
    Check that username and password are defined

    Arguments:
        username: Username
        password: Password

    Raises:
        MissingCredentialsError: If either `username` of `password` is `None`.
    """

    if token_path is not None and (
        username is not None or password is not None
    ):
        raise ValueError("Cannot set both token and username/password.")

    if not token_path and not username:
        message = (
            "FRACTAL_USER variable not defined."
            "\nPossible options: \n"
            + "    1. Set --user argument;\n"
            + "    2. Define FRACTAL_USER in a .fractal.env file;\n"
            + "    3. Define FRACTAL_USER as an environment variable."
        )
        raise MissingCredentialsError(message)
    if not token_path and not password:
        message = (
            "FRACTAL_PASSWORD variable not defined."
            "\nPossible options: \n"
            + "    1. Set --password argument;\n"
            + "    2. Define FRACTAL_PASSWORD in a .fractal.env file;\n"
            + "    3. Define FRACTAL_PASSWORD as an environment variable."
        )
        raise MissingCredentialsError(message)
    return (username, password)


def handle(cli_args: list[str] = argv):
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
        # no command provided. Print help and exit 1
        parser_main.print_help()
        exit(1)

    try:
        # Make a copy of vars(args), and remove cmd (which is not a relevant
        # argument for functions called with **kwargs)
        kwargs = vars(args).copy()
        kwargs.pop("cmd")
        fractal_server = (
            kwargs.pop("fractal_server") or settings.FRACTAL_SERVER
        )

        if args.cmd == "version":
            interface = handler(fractal_server, **kwargs)
        else:
            # Extract (and remove) username/password for AuthClient from kwargs
            username = kwargs.pop("user") or settings.FRACTAL_USER
            password = kwargs.pop("password") or settings.FRACTAL_PASSWORD
            token_path = (
                kwargs.pop("token_path") or settings.FRACTAL_TOKEN_PATH
            )
            _validate_credentials(
                username=username,
                password=password,
                token_path=token_path,
            )
            if token_path is not None:
                with open(token_path) as f:
                    token = f.read().strip("\n")
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
    exit(interface.retcode)


if __name__ == "__main__":
    main()
