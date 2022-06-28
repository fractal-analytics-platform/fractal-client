from os import getenv

from dotenv import load_dotenv
from pydantic import BaseSettings


def fail_getenv(key):
    value = getenv(key, None)
    if value is None:
        raise ValueError(f"Must provide environment variable {key}")
    return value


load_dotenv(".fractal.env")


__VERSION__ = "0.1.0"


class Settings(BaseSettings):
    PROJECT_NAME: str = "Fractal client"
    PROJECT_VERSION: str = __VERSION__

    FRACTAL_USER: str = fail_getenv("FRACTAL_USER")
    FRACTAL_PASSWORD: str = fail_getenv("FRACTAL_PASSWORD")

    FRACTAL_SERVER: str = getenv("FRACTAL_SERVER", "http://localhost:8000")

    BASE_URL: str = f"{FRACTAL_SERVER}/api/v1"


settings = Settings()
