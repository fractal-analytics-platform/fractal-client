import logging
import sys
from datetime import datetime
from pathlib import Path

import jwt

MIN_TOKEN_TTL = 10.0


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


def read_and_refresh_token(path: str) -> str:
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
