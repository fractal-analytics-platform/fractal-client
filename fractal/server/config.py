from enum import Enum
from os import getenv

from dotenv import load_dotenv
from pydantic import BaseSettings


def fail_getenv(key):
    value = getenv(key, None)
    if value is None:
        raise ValueError(f"Must provide environment variable {key}")
    return value


load_dotenv(".fractal.env")


class DeploymentType(str, Enum):
    production = "production"
    staging = "staging"
    testing = "testing"
    development = "development"


class Settings(BaseSettings):
    PROJECT_NAME: str = "Fractal Server"
    PROJECT_VERSION: str = "0.1.0"
    DEPLOYMENT_TYPE: DeploymentType = DeploymentType(
        fail_getenv("DEPLOYMENT_TYPE")
    )

    # AUTH
    LDAP_SERVER: str | None = getenv("LDAP_SERVER", None)
    LDAP_SSL: bool = getenv("LDAP_SSL", "1") == "1"
    JWT_EXPIRE_MINUTES: int = int(getenv("JWT_EXPIRE_MINUTES", default=3))
    JWT_SECRET_KEY: str = fail_getenv("JWT_SECRET_KEY")
    JWT_ALGORITHM: str = "HS256"


settings = Settings()
