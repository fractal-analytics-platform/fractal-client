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
from enum import Enum
from os import environ
from os import getenv
from os.path import abspath
from typing import List
from typing import Optional

from dotenv import load_dotenv
from pydantic import BaseModel
from pydantic import BaseSettings
from pydantic import Field
from pydantic import root_validator


def fail_getenv(key):
    value = getenv(key, None)
    if value is None:
        raise ValueError(f"Must provide environment variable {key}")
    return value


load_dotenv(".fractal_server.env")


class OAuthClient(BaseModel):
    CLIENT_NAME: str
    CLIENT_ID: str
    CLIENT_SECRET: str

    AUTHORIZE_ENDPOINT: Optional[str]
    ACCESS_TOKEN_ENDPOINT: Optional[str]
    REFRESH_TOKEN_ENDPOINT: Optional[str]
    REVOKE_TOKEN_ENDPOINT: Optional[str]


__VERSION__ = "0.1.0"


class DeploymentType(str, Enum):
    PRODUCTION = "production"
    STAGING = "staging"
    TESTING = "testing"
    DEVELOPMENT = "development"


class Settings(BaseSettings):
    PROJECT_NAME: str = "Fractal Server"
    PROJECT_VERSION: str = __VERSION__
    DEPLOYMENT_TYPE: DeploymentType = DeploymentType(
        fail_getenv("DEPLOYMENT_TYPE")
    )

    ###########################################################################
    # AUTH
    ###########################################################################

    OAUTH_CLIENTS: List[OAuthClient] = Field(default_factory=list)

    # JWT TOKEN
    JWT_EXPIRE_SECONDS: int = int(getenv("JWT_EXPIRE_SECONDS", default=180))
    JWT_SECRET_KEY: str = fail_getenv("JWT_SECRET_KEY")

    ###########################################################################
    # DATABASE
    ###########################################################################
    DB_ENGINE: str = getenv("DB_ENGINE", "sqlite")

    if DB_ENGINE == "postgres":
        POSTGRES_USER: str = getenv("POSTGRES_USER", "root")
        POSTGRES_PASSWORD = getenv("POSTGRES_PASSWORD", "password")
        POSTGRES_SERVER: str = getenv("POSTGRES_SERVER", "localhost")
        POSTGRES_PORT: str = getenv("POSTGRES_PORT", "5432")
        POSTGRES_DB: str = getenv("POSTGRES_DB", "test_db")

        DATABASE_URL = (
            f"postgresql+asyncpg://{POSTGRES_USER}:{POSTGRES_PASSWORD}"
            f"@{POSTGRES_SERVER}:{POSTGRES_PORT}/{POSTGRES_DB}"
        )
    elif DB_ENGINE == "sqlite":
        SQLITE_PATH: str = getenv("SQLITE_PATH", "")

        DATABASE_URL = (
            "sqlite+aiosqlite:///"
            f"{abspath(SQLITE_PATH) if SQLITE_PATH else SQLITE_PATH}"
        )

    @property
    def DB_ECHO(self):
        db_echo = bool(
            int(
                getenv(
                    "DB_ECHO",
                    self.DEPLOYMENT_TYPE != DeploymentType.PRODUCTION,
                )
            )
        )
        return db_echo

    ###########################################################################
    # FRACTAL SPECIFIC
    ###########################################################################
    DATA_DIR_ROOT: str = fail_getenv("DATA_DIR_ROOT")

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


settings = Settings()
