from pathlib import Path

import jwt
from httpx import AsyncClient
from jwt.exceptions import ExpiredSignatureError

from .config import settings


class AuthenticationError(ValueError):
    pass


class AuthToken:
    def __init__(
        self,
        client: AsyncClient,
        username: str,
        password: str,
        slurm_user: str,
    ):
        self.client = client
        self.username = username
        self.password = password
        self.slurm_user = slurm_user

        try:
            with open(f"{settings.FRACTAL_CACHE_PATH}/session", "r") as f:
                self.token = f.read()
        except FileNotFoundError:
            pass

    async def _get_fresh_token(self):
        data = dict(
            username=self.username,
            password=self.password,
            slurm_user=self.slurm_user,
        )
        res = await self.client.post(
            f"{settings.FRACTAL_SERVER}/auth/token/login", data=data
        )
        if res.status_code != 200:
            raise AuthenticationError(
                "Error: could not obtain token. Is the user registered?\n"
                f"{res.json()}\n"
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

    async def header(self):
        token = await self.__call__()
        return dict(Authorization=f"Bearer {token}")

    async def __call__(self):
        if self.expired:
            await self._get_fresh_token()
        return self.token


class AuthClient:
    def __init__(
        self,
        username: str,
        password: str,
        slurm_user: str,
    ):
        self.auth = None
        self.client = None
        self.username = username
        self.password = password
        self.slurm_user = slurm_user

    async def __aenter__(self):
        self.client = AsyncClient()
        self.auth = AuthToken(
            client=self.client,
            username=self.username,
            password=self.password,
            slurm_user=self.slurm_user,
        )
        return self

    async def __aexit__(self, exc_type, exc_value, traceback):
        await self.client.aclose()

    async def get(self, *args, **kwargs):
        return await self.client.get(
            headers=await self.auth.header(), *args, **kwargs
        )

    async def post(self, *args, **kwargs):
        return await self.client.post(
            headers=await self.auth.header(), *args, **kwargs
        )

    async def patch(self, *args, **kwargs):
        return await self.client.patch(
            headers=await self.auth.header(), *args, **kwargs
        )

    async def delete(self, *args, **kwargs):
        return await self.client.delete(
            headers=await self.auth.header(), *args, **kwargs
        )
