from pathlib import Path

import jwt
from httpx import AsyncClient
from jwt.exceptions import ExpiredSignatureError

from .config import settings


class AuthToken:
    def __init__(self, client: AsyncClient):
        self.client = client

        try:
            with open(settings.SESSION_CACHE_PATH, "r") as f:
                self.token = f.read()
        except FileNotFoundError:
            pass

    async def _get_fresh_token(self):
        data = dict(
            username=settings.FRACTAL_USER,
            password=settings.FRACTAL_PASSWORD,
        )
        res = await self.client.post(
            f"{settings.FRACTAL_SERVER}/auth/token/login", data=data
        )
        if res.status_code != 200:
            raise ValueError("ERROR! (probably the user was not registered)")
        raw_token = res.json()
        self.token = raw_token["access_token"]
        with open(Path(settings.SESSION_CACHE_PATH).expanduser(), "w") as f:
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

    async def header(self):
        token = await self.__call__()
        return dict(Authorization=f"Bearer {token}")

    async def __call__(self):
        if self.expired:
            await self._get_fresh_token()
        return self.token
