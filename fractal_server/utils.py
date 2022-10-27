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
from warnings import warn as _warn

from .config import Settings
from .syringe import Inject


def warn(message, settings: Settings = Inject(Settings)):
    """
    Make sure that warnings do not make their way to staing and production
    """
    if settings.DEPLOYMENT_TYPE in ["testing", "development"]:
        _warn(message, RuntimeWarning)
    else:
        raise RuntimeError(message)


def slugify(value: str):
    return value.lower().replace(" ", "_")
