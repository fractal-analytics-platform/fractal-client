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
from os import getenv
from pathlib import Path

from dotenv import load_dotenv

load_dotenv(".fractal.env")


class Settings:
    def __init__(self):

        self.FRACTAL_USER: str | None = getenv("FRACTAL_USER")
        self.FRACTAL_PASSWORD: str | None = getenv("FRACTAL_PASSWORD")
        self.FRACTAL_TOKEN_PATH: str | None = getenv("FRACTAL_TOKEN")

        self.FRACTAL_SERVER: str = getenv("FRACTAL_SERVER")
        self.FRACTAL_CACHE_PATH: str = getenv(
            "FRACTAL_CACHE_PATH", str(Path.home() / ".cache/fractal")
        )


settings = Settings()
