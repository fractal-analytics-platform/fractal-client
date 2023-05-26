from json import JSONDecodeError
from pathlib import Path

import jwt
from httpx import Client
from jwt.exceptions import ExpiredSignatureError

from .config import settings


class AuthenticationError(ValueError):
    pass


class AuthToken:
    def __init__(
        self,
        client: Client,
        username: str,
        password: str,
    ):
        self.client = client
        self.username = username
        self.password = password

        try:
            with open(f"{settings.FRACTAL_CACHE_PATH}/session", "r") as f:
                self.token = f.read()
        except FileNotFoundError:
            pass

    def _get_fresh_token(self):
        data = dict(
            username=self.username,
            password=self.password,
        )
        res = self.client.post(
            f"{settings.FRACTAL_SERVER}/auth/token/login", data=data
        )
        if res.status_code != 200:
            try:
                data = res.json()
            except JSONDecodeError:
                data = ""
            raise AuthenticationError(
                "Error: could not obtain token. Is the user registered?\n"
                f"{data}\n"
            )
        raw_token = res.json()
        self.token = raw_token["access_token"]

        # Create cache folder, if needed
        cache_dir = Path(f"{settings.FRACTAL_CACHE_PATH}").expanduser()
        cache_dir.mkdir(parents=True, exist_ok=True)

        # Write token in cache_file
        cache_file = cache_dir / "session"
        with cache_file.open("w") as f:
            f.write(self.token)

    @property
    def expired(self):
        try:
            jwt.decode(
                jwt=self.token,
                requires=["exp"],
                options={
                    "verify_signature": False,
                    "verify_exp": True,
                },
            )
            return False
        except AttributeError:
            return True
        except ExpiredSignatureError:
            return True

    def header(self):
        token = self.__call__()
        return dict(Authorization=f"Bearer {token}")

    def __call__(self):
        if self.expired:
            self._get_fresh_token()
        return self.token


class AuthClient:
    def __init__(
        self,
        username: str,
        password: str,
    ):
        self.auth = None
        self.client = None
        self.username = username
        self.password = password

    def __enter__(self):
        self.client = Client()
        self.auth = AuthToken(
            client=self.client,
            username=self.username,
            password=self.password,
        )
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.client.close()

    def get(self, *args, **kwargs):
        return self.client.get(headers=self.auth.header(), *args, **kwargs)

    def post(self, *args, **kwargs):
        return self.client.post(headers=self.auth.header(), *args, **kwargs)

    def patch(self, *args, **kwargs):
        return self.client.patch(headers=self.auth.header(), *args, **kwargs)

    def delete(self, *args, **kwargs):
        return self.client.delete(headers=self.auth.header(), *args, **kwargs)
