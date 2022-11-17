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
from os import getenv

from dotenv import load_dotenv
from pydantic import BaseSettings

from . import __VERSION__


def fail_getenv(key):
    value = getenv(key, None)
    if value is None:
        raise ValueError(f"Must provide environment variable {key}")
    return value


load_dotenv(".fractal.env")


class Settings(BaseSettings):
    PROJECT_NAME: str = "Fractal client"
    PROJECT_VERSION: str = __VERSION__

    FRACTAL_USER: str = fail_getenv("FRACTAL_USER")
    FRACTAL_PASSWORD: str = fail_getenv("FRACTAL_PASSWORD")
    SLURM_USER: str = fail_getenv("SLURM_USER")

    FRACTAL_SERVER: str = getenv("FRACTAL_SERVER", "http://localhost:8000")

    BASE_URL: str = f"{FRACTAL_SERVER}/api/v1"
    FRACTAL_CACHE_PATH: str = "~/.cache/fractal"


settings = Settings()
