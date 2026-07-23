import logging
import sys
from datetime import datetime
from json import JSONDecodeError
from pathlib import Path
from typing import Any
from typing import Self

import jwt
from httpx2 import Client
from httpx2 import Response

from fractal_client.config import settings

from ._info import AuthInfo

logging.getLogger("httpx2").setLevel(logging.WARNING)

MIN_TOKEN_TTL = 10.0


def _debug_request(verb: str, url):
    log = f"Sending HTTP request {verb} {url}"
    logging.debug(log)


class AuthenticationError(ValueError):
    pass


def _get_validity_seconds(token: str) -> int | float:
    """
    Token validity time left (in seconds).

    Returns a negative number when the token cannot be decoded or does
    not contain the "exp" claim.

    Note that `HS256` is the default algorithm in `fastapi-users`, see
    https://github.com/fastapi-users/fastapi-users/blob/master/fastapi_users/jwt.py
    """
    try:
        claims = jwt.decode(
            token,
            algorithms=["HS256"],
            options={"verify_signature": False},
        )
        expiration = datetime.fromtimestamp(claims["exp"])
        now = datetime.now()
        validity_seconds = (expiration - now).total_seconds()
        logging.debug(
            (
                f"Token validity time {validity_seconds:.1f} "
                f"seconds (expiration: {expiration})"
            )
        )
        return validity_seconds
    except Exception as e:
        logging.debug(f"An exception took place: {str(e)}")
        return -1


def _is_token_valid(token: str) -> bool:
    return _get_validity_seconds(token) > MIN_TOKEN_TTL


def _read_and_refresh_token(path: str) -> str:
    if Path(path).exists():
        with open(path) as f:
            token: str = f.read().strip()
        source_info = f"at {path}"
        update_token = False
    else:
        prompt_msg = (
            f"No token was found at {path}, and no other authentication method "
            "is available.\n"
            f"Paste a valid token here and it will be written to {path}: "
        )
        token: str = input(prompt_msg)
        source_info = "provided"
        update_token = True

    if not _is_token_valid(token):
        prompt_msg = (
            f"\nThe token {source_info} is invalid or expired/expiring.\n"
            "Please get a fresh token and paste it here: "
        )
        token = input(prompt_msg)
        if not _is_token_valid(token):
            sys.exit("\nThe token provided is invalid or expired/expiring. Exit.")
        update_token = True

    if update_token:
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w") as f:
            f.write(token)

    return token


class AuthClient:
    fractal_server: str
    auth_info: AuthInfo
    client: Client
    token: str

    def __init__(self: Self, *, fractal_server: str, auth_info: AuthInfo):
        self.fractal_server = fractal_server
        self.auth_info: AuthInfo = auth_info
        self.client: Client = Client()

    def __enter__(self: Self):
        if self.auth_info.use_basic_auth:
            self.token: str = self.get_token_from_backend()
        elif self.auth_info.token_path:
            self.token: str = _read_and_refresh_token(self.auth_info.token_path)
        else:
            self.token: str = _read_and_refresh_token(settings.default_token_path)

        return self

    def __exit__(self: Self, *args):
        self.client.close()

    def get_token_from_backend(self: Self) -> str:
        res: Response = self.client.post(
            f"{self.fractal_server}/auth/token/login/",
            data=dict(
                username=self.auth_info.user,
                password=self.auth_info.password,
            ),
        )
        if res.status_code != 200:
            data: str = res.text
            raise AuthenticationError(
                f"Error at {res.request.url}.\n"
                f"Status code: {res.status_code}.\n"
                f"Response data: {data}.\n"
            )
        try:
            raw_token: dict[str, Any] = res.json()
        except JSONDecodeError:
            print("Error while parsing")  # FIXME
        return raw_token["access_token"]

    @property
    def auth_headers(self: Self) -> dict[str, str]:
        return {"Authorization": f"Bearer {self.token}"}

    def _get_url(self: Self, relative_url: str) -> str:
        relative_url_no_leading_slash = relative_url.lstrip("/")
        return f"{self.fractal_server}/{relative_url_no_leading_slash}"

    def get(self: Self, relative_url: str):
        url = self._get_url(relative_url)
        _debug_request("GET", url)
        return self.client.get(url=url, headers=self.auth_headers)

    def post(self: Self, relative_url: str, **kwargs):
        url = self._get_url(relative_url)
        _debug_request("POST", url)
        return self.client.post(url=url, headers=self.auth_headers, **kwargs)

    def patch(self: Self, relative_url: str, **kwargs):
        url = self._get_url(relative_url)
        _debug_request("PATCH", url)
        return self.client.patch(url=url, headers=self.auth_headers, **kwargs)

    def delete(self: Self, relative_url: str):
        url = self._get_url(relative_url)
        _debug_request("DELETE", url)
        return self.client.delete(url=url, headers=self.auth_headers)
