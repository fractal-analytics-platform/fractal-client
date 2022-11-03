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
import asyncio
import logging
from sys import argv
from typing import List

from httpx import AsyncClient
from httpx import ConnectError

from . import cmd
from .authclient import AuthClient
from .authclient import AuthenticationError
from .config import settings
from .interface import PrintInterface
from .parser import parser_main


logging.basicConfig(level=logging.DEBUG)


async def handle(cli_args: List[str] = argv):
    args = parser_main.parse_args(cli_args[1:])
    logging.debug(args)

    if args.cmd:
        handler = getattr(cmd, args.cmd)
    else:
        # no command provided. Print help and exit 1
        parser_main.print_help()
        exit(1)

    try:
        if args.cmd in ["version", "register"]:
            async with AsyncClient() as client:
                interface = await handler(client, **vars(args))
        else:
            async with AuthClient(
                username=args.user or settings.FRACTAL_USER,
                password=args.password or settings.FRACTAL_PASSWORD,
                slurm_user=args.slurm_user or settings.SLURM_USER,
            ) as client:
                interface = await handler(client, **vars(args))
    except AuthenticationError as e:
        return PrintInterface(retcode=1, data=e.args[0])
    except ConnectError as e:
        return PrintInterface(retcode=1, data=e.args[0])

    return interface


async def main():
    interface = await handle()
    interface.show()
    exit(interface.retcode)


if __name__ == "__main__":
    asyncio.run(main())
