import logging
import sys
from datetime import datetime
from json import JSONDecodeError
from pathlib import Path
from typing import Self

import jwt
from httpx2 import Client

from fractal_client.config import settings

from .info import AuthInfo

logging.getLogger("httpx2").setLevel(logging.WARNING)

MIN_TOKEN_TTL = 10.0


def debug_request(verb: str, url):
    log = f"Sending HTTP request {verb} {url}"
    logging.debug(log)


class AuthenticationError(ValueError):
    pass


def get_ttl(token: str) -> int | float:
    """
    Token validity time left (in seconds).

    Returns a negative number when the token cannot be decoded or does
    not contain the "exp" claim.
    """
    try:
        claims = jwt.decode(
            token,
            algorithms=["HS256"],
            options={"verify_signature": False},
        )
        expiration = datetime.fromtimestamp(claims["exp"])
        now = datetime.now()
        ttl = expiration - now
        logging.debug(
            (
                f"Token validity time {ttl.total_seconds():.1f} "
                f"seconds (expiration: {expiration})"
            )
        )
        return ttl.total_seconds()
    except Exception as e:
        logging.debug(f"An exception took place: {str(e)}")
        return -1


def is_token_valid(token: str) -> bool:
    return get_ttl(token) > MIN_TOKEN_TTL


def read_custom_token(path: str) -> str:
    if not Path(path).exists():
        raise  # FIXME
    with open(path) as f:
        token = f.read().strip()

    if not is_token_valid(token):
        prompt_msg = (
            f"\nThe token stored at {path} is invalid or expired/expiring.\n"
            "Please get a fresh token and paste it here: "
        )
        token = input(prompt_msg)
        if not is_token_valid(token):
            sys.exit("\nThe provided token is invalid or expired/expiring. Exit.")

    return token


def read_and_write_standard_token() -> str:
    path: str = settings.default_token_path
    if not Path(path).exists():
        sys.exit(
            f"No token found at {path}, and no other "
            "authentication method is available."
        )
    with open(path) as f:
        token = f.read().strip()

    if not is_token_valid(token):
        prompt_msg = (
            f"\nThe token stored at {path} is invalid or expired/expiring.\n"
            "Please get a fresh token and paste it here: "
        )
        token = input(prompt_msg)
        if not is_token_valid(token):
            sys.exit("\nThe provided token is invalid or expired/expiring. Exit.")
    with open(path, "w") as f:
        f.write(token)

    return token


class AuthClient:
    fractal_server: str
    auth_info: AuthInfo
    client: Client
    token: str

    def __init__(self, *, fractal_server: str, auth_info: AuthInfo):
        self.fractal_server = fractal_server
        self.auth_info: AuthInfo = auth_info
        self.client = Client()

    def __enter__(self: Self):
        if self.auth_info.use_basic_auth:
            self.token = self.get_token_from_backend()
        elif self.auth_info.token_path:
            self.token: str = read_custom_token(self.auth_info.token_path)
        else:
            self.token = read_and_write_standard_token()

        return self

    def __exit__(self, *args):
        self.client.close()

    def get_token_from_backend(self) -> str:
        res = self.client.post(
            f"{self.fractal_server}/auth/token/login/",
            data=dict(
                username=self.auth_info.user,
                password=self.auth_info.password,
            ),
        )
        if res.status_code != 200:
            data = res.text
            raise AuthenticationError(
                f"Error at {res.request.url}.\n"
                f"Status code: {res.status_code}.\n"
                f"Response data: {data}.\n"
            )
        try:
            raw_token = res.json()
        except JSONDecodeError:
            raise  # FIXME
        return raw_token["access_token"]

    @property
    def auth_headers(self) -> dict[str, str]:
        return {"Authorization": f"Bearer {self.token}"}

    def _get_url(self, relative_url: str) -> str:
        relative_url_no_leading_slash = relative_url.lstrip("/")
        return f"{self.fractal_server}/{relative_url_no_leading_slash}"

    def get(self, relative_url: str):
        url = self._get_url(relative_url)
        debug_request("GET", url)
        return self.client.get(url=url, headers=self.auth_headers)

    def post(self, relative_url: str, **kwargs):
        url = self._get_url(relative_url)
        debug_request("POST", url)
        return self.client.post(url=url, headers=self.auth_headers, **kwargs)

    def patch(self, relative_url: str, **kwargs):
        url = self._get_url(relative_url)
        debug_request("PATCH", url)
        return self.client.patch(url=url, headers=self.auth_headers, **kwargs)

    def delete(self, relative_url: str):
        url = self._get_url(relative_url)
        debug_request("DELETE", url)
        return self.client.delete(url=url, headers=self.auth_headers)
