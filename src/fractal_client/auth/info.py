from typing import Self

from pydantic import BaseModel


class AuthInfo(BaseModel):
    user: str | None
    password: str | None
    token_path: str | None

    @property
    def use_basic_auth(self: Self) -> bool:
        return bool(self.user) and bool(self.password)
