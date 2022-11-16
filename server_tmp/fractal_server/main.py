"""
Copyright 2022 (C) Friedrich Miescher Institute for Biomedical Research and
University of Zurich

Original authors:
Jacopo Nespolo <jacopo.nespolo@exact-lab.it>

This file is part of Fractal and was originally developed by eXact lab S.r.l.
<exact-lab.it> under contract with Liberali Lab from the Friedrich Miescher
Institute for Biomedical Research and Pelkmans Lab from the University of
Zurich.
"""
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from .config import get_settings
from .syringe import Inject


def collect_routers(app: FastAPI) -> None:
    from .app.api import router_default
    from .app.api import router_v1
    from .app.security import auth_router

    app.include_router(router_default, prefix="/api")
    app.include_router(router_v1, prefix="/api/v1")
    app.include_router(auth_router, prefix="/auth", tags=["auth"])


def check_settings():
    """
    Check and register the settings
    """
    settings = Inject(get_settings)
    settings.check()


async def __on_startup():
    check_settings()


def start_application() -> FastAPI:
    app = FastAPI()
    collect_routers(app)

    app.add_middleware(
        CORSMiddleware,
        allow_origins=["http://127.0.0.1:5173", "http://localhost:5173"],
        allow_methods=["GET", "POST", "PUT", "PATCH", "DELETE", "OPTIONS"],
        allow_headers=[
            "set-cookie",
            "Set-Cookie",
            "Content-Type",
            "Access-Control-Allow-Headers",
            "X-Requested-With",
        ],
        allow_credentials=True,
    )

    return app


app = start_application()


@app.on_event("startup")
async def on_startup():
    await __on_startup()
