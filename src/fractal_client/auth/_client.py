import logging
from json import JSONDecodeError
from typing import Any
from typing import Self

from httpx2 import Client
from httpx2 import Response

from fractal_client.auth._args import AuthInfo
from fractal_client.auth._token_utils import read_and_refresh_token
from fractal_client.config import settings

logging.getLogger("httpx2").setLevel(logging.WARNING)


def _debug_request(verb: str, url):
    log = f"Sending HTTP request {verb} {url}"
    logging.debug(log)


class AuthenticationError(ValueError):
    pass


class AuthClient:
    fractal_server: str
    auth_info: AuthInfo
    client: Client
    token: str

    def __init__(self: Self, *, fractal_server: str, auth_info: AuthInfo):
        self.fractal_server: str = fractal_server
        self.auth_info: AuthInfo = auth_info
        self.client: Client = Client()

    def __enter__(self: Self):
        if self.auth_info.use_basic_auth:
            self.token: str = self._get_token_from_backend()
        elif self.auth_info.token_path:
            self.token: str = read_and_refresh_token(self.auth_info.token_path)
        else:
            self.token: str = read_and_refresh_token(settings.default_token_path)

        return self

    def __exit__(self: Self, *args):
        self.client.close()

    def _get_token_from_response(self: Self, response: Response) -> str:
        try:
            response_data: dict[str, Any] = response.json()
            token_value: str = response_data["access_token"]
            return token_value
        except (JSONDecodeError, KeyError) as e:
            raise AuthenticationError(
                "Could not extract token from Fractal-backend response "
                f"({type(e).__name__})."
            ) from e

    def _get_token_from_backend(self: Self) -> str:
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
        return self._get_token_from_response(res)

    @property
    def _auth_headers(self: Self) -> dict[str, str]:
        return {"Authorization": f"Bearer {self.token}"}

    def _get_url(self: Self, relative_url: str) -> str:
        relative_url_no_leading_slash = relative_url.lstrip("/")
        return f"{self.fractal_server}/{relative_url_no_leading_slash}"

    def get(self: Self, relative_url: str):
        url = self._get_url(relative_url)
        _debug_request("GET", url)
        return self.client.get(url=url, headers=self._auth_headers)

    def post(self: Self, relative_url: str, **kwargs):
        url = self._get_url(relative_url)
        _debug_request("POST", url)
        return self.client.post(url=url, headers=self._auth_headers, **kwargs)

    def patch(self: Self, relative_url: str, **kwargs):
        url = self._get_url(relative_url)
        _debug_request("PATCH", url)
        return self.client.patch(url=url, headers=self._auth_headers, **kwargs)

    def delete(self: Self, relative_url: str):
        url = self._get_url(relative_url)
        _debug_request("DELETE", url)
        return self.client.delete(url=url, headers=self._auth_headers)
