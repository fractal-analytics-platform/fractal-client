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

from pydantic import BaseSettings


class Settings(BaseSettings):
    RUNNER_CONFIG: str = getenv("RUNNER_CONFIG", "local")
    RUNNER_LOG_DIR: str = getenv("RUNNER_LOG_DIR", "logs")
    RUNNER_DEFAULT_EXECUTOR: str = getenv("RUNNER_DEFAULT_EXECUTOR", "cpu-low")
    RUNNER_MONITORING: bool = int(getenv("RUNNER_MONITORING", 1))


settings = Settings()
