import sys
from pathlib import Path

from httpx2 import Client
from httpx2._models import Response

from fractal_client.auth._token_utils import _get_token_hint
from fractal_client.auth._token_utils import _is_token_valid
from fractal_client.config import settings
from fractal_client.interface import Interface


def auth_check_token(*, fractal_server: str, token_path: str | None) -> Interface:
    path = Path(token_path or settings.default_token_path)
    if not path.exists():
        sys.exit(f"File not found at {path}.")
    token = Path(path).read_text().strip()
    try:
        with Client() as client:
            url = f"{fractal_server}/auth/current-user/"
            res: Response = client.get(
                url, headers={"Authorization": f"Bearer {token}"}
            )
            if res.status_code != 200:
                raise ValueError(f"Server responded with {res}")
            user = res.json()
            user_email = user["email"]
        return Interface(
            retcode=0,
            data=f"Valid token for {user_email} found at {path}.",
        )
    except Exception as e:
        return Interface(
            retcode=1,
            data=f"Token verification failed. Original error: {str(e)}",
        )


def auth_set_token(*, token_path: str | None) -> Interface:
    path = Path(token_path or settings.default_token_path)
    hint = _get_token_hint()
    prompt_msg = f"Paste a valid token here{hint}, and it will be written to {path}: "
    token = input(prompt_msg).strip()
    if _is_token_valid(token):
        path.parent.mkdir(exist_ok=True, parents=True)
        path.write_text(token)
        return Interface(
            retcode=0,
            data=f"Token written to {path}.",
        )
    else:
        return Interface(
            retcode=1,
            data="The token provided is invalid or expired/expiring. Exit.",
        )


def auth_clear_token(*, token_path: str | None) -> Interface:
    path = Path(token_path or settings.default_token_path)
    if path.exists():
        path.unlink()
        return Interface(
            retcode=0,
            data=f"Removed {path}.",
        )
    else:
        return Interface(
            retcode=1,
            data=f"File not found at {path}, exit.",
        )
