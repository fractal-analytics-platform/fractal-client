"""
Copyright 2022 (C) Friedrich Miescher Institute for Biomedical Research and
University of Zurich

Original authors:
Jacopo Nespolo <jacopo.nespolo@exact-lab.it>
Tommaso Comparin <tommaso.comparin@exact-lab.it>

This file is part of Fractal and was originally developed by eXact lab S.r.l.
<exact-lab.it> under contract with Liberali Lab from the Friedrich Miescher
Institute for Biomedical Research and Pelkmans Lab from the University of
Zurich.
"""
import logging
from os import getenv
from pathlib import Path
from typing import Optional

from dotenv import load_dotenv

from . import __VERSION__

load_dotenv(".fractal.env")


class Settings:
    def __init__(self):
        self.PROJECT_NAME: str = "Fractal client"
        self.PROJECT_VERSION: str = __VERSION__

        self.FRACTAL_LOGGING_LEVEL: int = logging.INFO

        self.FRACTAL_USER: Optional[str] = getenv("FRACTAL_USER")
        self.FRACTAL_PASSWORD: Optional[str] = getenv("FRACTAL_PASSWORD")

        self.FRACTAL_SERVER: str = getenv(
            "FRACTAL_SERVER", "http://localhost:8000"
        )

        self.BASE_URL: str = f"{self.FRACTAL_SERVER}/api/v1"
        self.FRACTAL_CACHE_PATH: str = str(Path.home() / ".cache/fractal")


settings = Settings()
