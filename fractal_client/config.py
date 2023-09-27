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

load_dotenv(".fractal.env")


class Settings:
    def __init__(self):

        self.FRACTAL_LOGGING_LEVEL: int = getenv(
            "FRACTAL_LOGGING_LEVEL", logging.INFO
        )

        self.FRACTAL_USER: Optional[str] = getenv("FRACTAL_USER")
        self.FRACTAL_PASSWORD: Optional[str] = getenv("FRACTAL_PASSWORD")

        self.FRACTAL_SERVER: str = getenv(
            "FRACTAL_SERVER", "http://localhost:8000"
        )
        self.FRACTAL_CACHE_PATH: str = getenv(
            "FRACTAL_CACHE_PATH", str(Path.home() / ".cache/fractal")
        )

    @property
    def BASE_URL(self) -> str:
        return f"{self.FRACTAL_SERVER}/api/v1"


settings = Settings()
