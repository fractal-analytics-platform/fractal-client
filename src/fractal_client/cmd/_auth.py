import sys
from pathlib import Path
from typing import Any

from httpx2 import Client
from httpx2._models import Response

from fractal_client.auth._token_utils import _get_token_hint
from fractal_client.auth._token_utils import _is_token_valid
from fractal_client.config import settings
from fractal_client.interface import Interface
from fractal_client.response import check_response


def auth_check_token(*, fractal_server: str, token_path: str | None) -> Interface:
    path = Path(token_path or settings.default_token_path)
    if not path.exists():
        sys.exit(f"File not found at {path}")
    token = Path(path).read_text().strip()
    with Client() as client:
        url = f"{fractal_server}/auth/current-user/"
        res: Response = client.get(url, headers={"Authorization": f"Bearer {token}"})
        user: dict[str, Any] = check_response(res, expected_status_code=200)
        user_email = user["email"]

    return Interface(retcode=0, data=f"Valid token for {user_email} found at {path}")


def auth_set_token(*, token_path: str | None) -> Interface:
    path = Path(token_path or settings.default_token_path)
    hint = _get_token_hint()
    prompt_msg = f"Paste a valid token here{hint}, and it will be written to {path}: "
    token = input(prompt_msg).strip()
    if _is_token_valid(token):
        path.parent.mkdir(exist_ok=True, parents=True)
        path.write_text(token)
    else:
        return Interface(
            retcode=1,
            data="The token provided is invalid or expired/expiring.",
        )

    return Interface(retcode=0, data="")


def auth_clear_token(*, token_path: str | None) -> Interface:
    path = Path(token_path or settings.default_token_path)
    if path.exists():
        path.unlink()
        return Interface(retcode=0, data=f"Removed {path}.")
    else:
        return Interface(retcode=1, data=f"File not found at {path}, exit.")
