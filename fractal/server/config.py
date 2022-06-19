from os import environ
from os import getenv

from dotenv import load_dotenv
from pydantic import BaseSettings


load_dotenv(".fractal.env")


class Settings(BaseSettings):
    PROJECT_NAME: str = "Fractal Server"
    PROJECT_VERSION: str = "0.1.0"

    # AUTH
    LDAP_SERVER: str | None = getenv("LDAP_SERVER", None)
    LDAP_SSL: bool = getenv("LDAP_SSL", "1") == "1"
    JWT_EXPIRE_MINUTES: int = int(getenv("JWT_EXPIRE_MINUTES", default=3))
    JWT_SECRET_KEY: str = environ["JWT_SECRET_KEY"]
    JWT_ALGORITHM: str = "HS256"


settings = Settings()
