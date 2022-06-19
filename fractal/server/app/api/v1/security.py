from fastapi import APIRouter
from fastapi import Depends
from fastapi import HTTPException
from fastapi import status
from fastapi.security import OAuth2PasswordRequestForm

from ...security import access_token_encode
from ...security import authenticate_user
from ...security import FailedAuthenticationException
from ...security import get_current_user
from ...security import Token
from ...security import User


router = APIRouter()


@router.post("/token", response_model=Token)
async def login_for_access_token(
    form_data: OAuth2PasswordRequestForm = Depends(),
):
    try:
        user = await authenticate_user(form_data.username, form_data.password)
    except FailedAuthenticationException:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Incorrect username or password",
        )

    access_token = access_token_encode(
        data={"sub": user.sub, "scopes": form_data.scopes},
    )
    return {"access_token": access_token, "token_type": "bearer"}


@router.get("/me", response_model=User)
async def who_am_i(current_user: User = Depends(get_current_user)):
    return current_user
