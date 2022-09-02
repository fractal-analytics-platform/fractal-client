"""
Copyright 2022 (C) Friedrich Miescher Institute for Biomedical Research and
University of Zurich

Original author(s):
Jacopo Nespolo <jacopo.nespolo@exact-lab.it>

This file is part of Fractal and was originally developed by eXact lab S.r.l.
<exact-lab.it> under contract with Liberali Lab from the Friedrich Miescher
Institute for Biomedical Research and Pelkmans Lab from the University of
Zurich.
"""
from fastapi import APIRouter
from fastapi import Depends

from ...db import AsyncSession
from ...db import get_db
from ...security import current_active_user
from ...security import User

router = APIRouter()


@router.post("/")
async def create_workflow(
    user: User = Depends(current_active_user),
    db: AsyncSession = Depends(get_db),
):
    """
    Add new global workflow
    """
    raise NotImplementedError


@router.post("/{wofklow_id}")
async def add_task(
    user: User = Depends(current_active_user),
    db: AsyncSession = Depends(get_db),
):
    """
    Add global task to global workflow
    """
    raise NotImplementedError


@router.get("/")
async def get_list_workflow(
    user: User = Depends(current_active_user),
    db: AsyncSession = Depends(get_db),
):
    """
    List public workflows
    """
    raise NotImplementedError
