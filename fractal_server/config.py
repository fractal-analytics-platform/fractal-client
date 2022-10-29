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
import logging
from os import environ
from os import getenv
from os.path import abspath
from pathlib import Path
from typing import List
from typing import Literal
from typing import Optional
from typing import TypeVar

from dotenv import load_dotenv
from pydantic import BaseModel
from pydantic import BaseSettings
from pydantic import Field
from pydantic import root_validator

import fractal_server


T = TypeVar("T")


load_dotenv(".fractal_server.env")


class OAuthClient(BaseModel):
    CLIENT_NAME: str
    CLIENT_ID: str
    CLIENT_SECRET: str

    AUTHORIZE_ENDPOINT: Optional[str]
    ACCESS_TOKEN_ENDPOINT: Optional[str]
    REFRESH_TOKEN_ENDPOINT: Optional[str]
    REVOKE_TOKEN_ENDPOINT: Optional[str]


class Settings(BaseSettings):
    class Config:
        case_sensitive = True

    PROJECT_NAME: str = "Fractal Server"
    PROJECT_VERSION: str = fractal_server.__VERSION__
    DEPLOYMENT_TYPE: Optional[
        Literal["production", "staging", "testing", "development"]
    ]

    ###########################################################################
    # AUTH
    ###########################################################################

    OAUTH_CLIENTS: List[OAuthClient] = Field(default_factory=list)

    # JWT TOKEN
    JWT_EXPIRE_SECONDS: int = 180
    JWT_SECRET_KEY: Optional[str]

    # COOKIE TOKEN
    COOKIE_EXPIRE_SECONDS: int = 86400

    @root_validator(pre=True)
    def collect_oauth_clients(cls, values):
        oauth_env_variable_keys = [
            key for key in environ.keys() if "OAUTH" in key
        ]
        clients_available = {
            var.split("_")[1] for var in oauth_env_variable_keys
        }

        values["OAUTH_CLIENTS"] = []
        for client in clients_available:
            prefix = f"OAUTH_{client}"
            oauth_client = OAuthClient(
                CLIENT_NAME=client,
                CLIENT_ID=getenv(f"{prefix}_CLIENT_ID", None),
                CLIENT_SECRET=getenv(f"{prefix}_CLIENT_SECRET", None),
                AUTHORIZE_ENDPOINT=getenv(
                    f"{prefix}_AUTHORIZE_ENDPOINT", None
                ),
                ACCESS_TOKEN_ENDPOINT=getenv(
                    f"{prefix}_ACCESS_TOKEN_ENDPOINT", None
                ),
                REFRESH_TOKEN_ENDPOINT=getenv(
                    f"{prefix}_REFRESH_TOKEN_ENDPOINT", None
                ),
                REVOKE_TOKEN_ENDPOINT=getenv(
                    f"{prefix}_REVOKE_TOKEN_ENDPOINT", None
                ),
            )
            values["OAUTH_CLIENTS"].append(oauth_client)
        return values

    ###########################################################################
    # DATABASE
    ###########################################################################
    DB_ENGINE: Literal["sqlite", "postgres"] = "sqlite"
    DB_ECHO: bool = False

    if DB_ENGINE == "postgres":
        POSTGRES_USER: str = Field()
        POSTGRES_PASSWORD: str = Field()
        POSTGRES_SERVER: str = "localhost"
        POSTGRES_PORT: str = "5432"
        POSTGRES_DB: str = Field()

    elif DB_ENGINE == "sqlite":
        SQLITE_PATH: Optional[str]

    @property
    def DATABASE_URL(self):
        if self.DB_ENGINE == "sqlite":
            if not self.SQLITE_PATH:
                raise ValueError("SQLITE_PATH path cannot be None")
            sqlite_path = (
                abspath(self.SQLITE_PATH)
                if self.SQLITE_PATH
                else self.SQLITE_PATH
            )
            return f"sqlite+aiosqlite:///{sqlite_path}"
        elif "postgres":
            pg_uri = (
                "postgresql+asyncpg://"
                f"{self.POSTGRES_USER}:{self.POSTGRES_PASSWORD}"
                f"@{self.POSTGRES_SERVER}:{self.POSTGRES_PORT}"
                f"/{self.POSTGRES_DB}"
            )
            return pg_uri

    @property
    def DATABASE_SYNC_URL(self):
        if self.DB_ENGINE == "sqlite":
            if not self.SQLITE_PATH:
                raise ValueError("SQLITE_PATH path cannot be None")
            return self.DATABASE_URL.replace("aiosqlite", "pysqlite")
        elif self.DB_ENGINE == "postgres":
            return self.DATABASE_URL.replace("asyncpg", "psycopg2")

    ###########################################################################
    # FRACTAL SPECIFIC
    ###########################################################################
    FRACTAL_ROOT: Optional[Path]  # Path = _DEVNULL
    RUNNER_BACKEND: str = "process"
    RUNNER_ROOT_DIR: Path = Path("artifacts")
    FRACTAL_LOGGING_LEVEL: int = logging.WARNING

    RUNNER_CONFIG: str = "local"
    RUNNER_DEFAULT_EXECUTOR: str = "cpu-low"

    # NOTE: we currently set RUNNER_MONITORING to False, due to
    # https://github.com/fractal-analytics-platform/fractal-server/issues/148
    # RUNNER_MONITORING: bool = int(getenv("RUNNER_MONITORING", 1))
    RUNNER_MONITORING: bool = False

    ###########################################################################
    # BUSINESS LOGIC
    ###########################################################################

    def check(self):
        """
        Make sure that mandatory variables are set

        This method must be called before the server starts
        """

        class StrictSettings(BaseSettings):
            class Config:
                extra = "allow"

            DEPLOYMENT_TYPE: Literal[
                "production", "staging", "testing", "development"
            ]
            JWT_SECRET_KEY: str
            DB_ENGINE: str = "sqlite"

            if DB_ENGINE == "postgres":
                POSTGRES_USER: str
                POSTGRES_PASSWORD: str
                POSTGRES_DB: str
            elif DB_ENGINE == "sqlite":
                SQLITE_PATH: str

            FRACTAL_ROOT: Path

        StrictSettings(**self.dict())


def get_settings(settings=Settings()) -> Settings:
    return settings
