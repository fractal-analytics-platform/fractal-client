from fastapi import APIRouter
from fastapi import Depends
from fastapi import HTTPException
from fastapi import Security
from fastapi import status
from fastapi.security import OAuth2PasswordRequestForm

from ...security import FailedAuthenticationException
from ...security import get_current_user
from ...security import Token
from ...security import token_authentication_login
from ...security import User


router = APIRouter()


@router.post("/token", response_model=Token)
async def login_for_access_token(
    form_data: OAuth2PasswordRequestForm = Depends(),
):
    try:
        return await token_authentication_login(
            form_data.username, form_data.password
        )
    except FailedAuthenticationException:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Incorrect username or password",
        )


@router.get("/me", response_model=User)
async def who_am_i(
    current_user: User = Security(get_current_user, scopes=["profile"])
):
    return current_user
