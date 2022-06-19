from fastapi import APIRouter
from fastapi import Depends
from fastapi import HTTPException
from fastapi import status
from fastapi.security import OAuth2PasswordRequestForm

from ...security import authenticate_user
from ...security import create_access_token
from ...security import FailedAuthenticationException
from ...security import Token


router = APIRouter()


@router.post("/token", response_model=Token)
async def login_for_access_token(
    form_data: OAuth2PasswordRequestForm = Depends(),
):
    try:
        user = await authenticate_user(form_data.username, form_data.password)
    except FailedAuthenticationException:
        raise HTTPException(
            status_code=status.HTTP401_UNAUTHORIZED,
            detail="Incorrect username or password",
        )

    access_token = create_access_token(
        data={"sub": user.sub, "scopes": form_data.scopes},
    )
    return {"access_token": access_token, "token_type": "bearer"}
