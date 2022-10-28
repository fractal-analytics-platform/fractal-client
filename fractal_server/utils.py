"""
Copyright 2022 (C) Friedrich Miescher Institute for Biomedical Research and
University of Zurich

Original authors:
Jacopo Nespolo <jacopo.nespolo@exact-lab.it>

This file is part of Fractal and was originally developed by eXact lab S.r.l.
<exact-lab.it> under contract with Liberali Lab from the Friedrich Miescher
Institute for Biomedical Research and Pelkmans Lab from the University of
Zurich.
"""
import asyncio
import logging
from collections.abc import MutableMapping
from datetime import datetime
from datetime import timezone
from pathlib import Path
from shlex import split as shlex_split
from typing import Any
from typing import Optional
from warnings import warn as _warn

from .config import Settings
from .syringe import Inject


def get_timestamp() -> datetime:
    return datetime.now(tz=timezone.utc)


def popget(d: MutableMapping, key: str, default: Any = None) -> Any:
    """
    Pop and return mapping item if possible or return default
    """
    try:
        return d.pop(key)
    except KeyError:
        return default


def warn(message):
    """
    Make sure that warnings do not make their way to staing and production
    """
    settings = Inject(Settings)
    if settings.DEPLOYMENT_TYPE in ["testing", "development"]:
        _warn(message, RuntimeWarning)
    else:
        raise RuntimeError(message)


def slugify(value: str):
    return value.lower().replace(" ", "_")


def set_logger(
    *,
    logger_name: str,
    level: int = logging.WARNING,
    log_file_path: Optional[Path] = None,
    formatter: Optional[logging.Formatter] = None,
) -> logging.Logger:
    """
    Set up and return a logger
    """
    logger = logging.getLogger(logger_name)
    logger.setLevel(level)
    if log_file_path:
        file_handler = logging.FileHandler(log_file_path, mode="a")
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    return logger


async def execute_command(
    *, cwd: Path, command: str, logger_name: Optional[str] = None
) -> str:
    """
    Execute arbitrary command

    If the command returns a return code different from zero, a RuntimeError
    containing the stderr is raised.

    Parameters
    ----------
    cwd : Path
        the working directory for the command execution
    command : str
        the command to execute

    Return
    ------
    stdout : str
        the stdout from the command execution
    """
    command_split = shlex_split(command)
    cmd, *args = command_split

    proc = await asyncio.create_subprocess_exec(
        cmd,
        *args,
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.PIPE,
        cwd=cwd,
    )
    stdout, stderr = await proc.communicate()
    logger = logging.getLogger(logger_name)
    logger.debug(f"Subprocess call to: {command}")
    if proc.returncode != 0:
        raise RuntimeError(stderr.decode("utf-8"))
    return stdout.decode("utf-8")
