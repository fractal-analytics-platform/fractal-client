import logging

from httpx import Client


logging.getLogger("httpx").setLevel(logging.WARNING)


def debug_request(verb: str, url: str, **kwargs):
    body = kwargs.get("json")
    log = f"\nFractal Client sending HTTP request to:\n    {verb} {url}"
    if body is not None:
        log += "\nRequest body:\n" + "\n".join(
            [f"    {k}: {v}" for k, v in body.items()]
        )
    logging.debug(log)


class AuthenticationError(ValueError):
    pass


class AuthClient:
    def __init__(
        self,
        *,
        fractal_server: str,
        username: str | None,
        password: str | None,
        token: str | None,
    ):
        self.fractal_server = fractal_server.rstrip("/")
        self.auth = None
        self.client = None
        self.username = username
        self.password = password
        self.token = token

    def __enter__(self):
        self.client = Client()
        if self.token is None:
            self.token = self._get_fresh_token()

        return self

    def __exit__(self, *args):
        self.client.close()

    def _get_fresh_token(self) -> str:
        res = self.client.post(
            f"{self.fractal_server}/auth/token/login/",
            data=dict(
                username=self.username,
                password=self.password,
            ),
        )
        if res.status_code != 200:
            data = res.text
            raise AuthenticationError(
                f"Error at {res.request.url}.\n"
                f"Status code: {res.status_code}.\n"
                f"Response data: {data}.\n"
            )
        raw_token = res.json()
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
        debug_request("POST", relative_url, **kwargs)
        return self.client.post(url=url, headers=self.auth_headers, **kwargs)

    def patch(self, relative_url: str, **kwargs):
        url = self._get_url(relative_url)
        debug_request("PATCH", relative_url, **kwargs)
        return self.client.patch(url=url, headers=self.auth_headers, **kwargs)

    def delete(self, relative_url: str):
        url = self._get_url(relative_url)
        debug_request("DELETE", url)
        return self.client.delete(url=url, headers=self.auth_headers)
