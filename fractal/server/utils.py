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
from collections.abc import MutableMapping
from typing import Any
from warnings import warn as _warn

from .config import DeploymentType
from .config import settings


def warn(message):
    """
    Make sure that warnings do not make their way to staing and production
    """
    if settings.DEPLOYMENT_TYPE in [
        DeploymentType.TESTING,
        DeploymentType.DEVELOPMENT,
    ]:
        _warn(message, RuntimeWarning)
    else:
        raise RuntimeError(message)


def popget(d: MutableMapping, key: str, default: Any = None) -> Any:
    """
    Pop and return mapping item if possible or return default
    """
    try:
        return d.pop(key)
    except KeyError:
        return default
