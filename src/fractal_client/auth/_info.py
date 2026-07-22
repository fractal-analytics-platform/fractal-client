import argparse as ap
import sys
from dataclasses import dataclass
from typing import Self

from fractal_client.config import settings

ERROR_MSG = (
    "Invalid authentication credentials. "
    "You should either set username&password or the token path.\n\n"
    "You can set these variables in multiple ways "
    "(see `fractal --help`):\n"
    "  1. Through command-line arguments (e.g. `-u`, `-p`, `--token-path`).\n"
    "  2. Through environment variables (e.g. `FRACTAL_USER`).\n"
    "  3. Through environment variables in a `.fractal.env` file.\n"
)


@dataclass(frozen=True)
class AuthInfo:
    user: str | None
    password: str | None
    token_path: str | None

    @property
    def use_basic_auth(self: Self) -> bool:
        return bool(self.user) and bool(self.password)


def get_auth_info(*, parser: ap.ArgumentParser, args: ap.Namespace) -> AuthInfo:
    user: str | None = args.user or settings.FRACTAL_SERVER
    password: str | None = args.password or settings.FRACTAL_PASSWORD
    token_path: str | None = args.token_path or settings.FRACTAL_TOKEN_PATH

    if bool(user) != bool(password):
        print(ERROR_MSG)
        sys.exit(1)

    if bool(user) and bool(token_path):
        print(ERROR_MSG)
        sys.exit(1)

    return AuthInfo(
        user=user,
        password=password,
        token_path=token_path,
    )
