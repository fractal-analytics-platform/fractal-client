"""
Copyright 2022 (C) Friedrich Miescher Institute for Biomedical Research and
University of Zurich

Original authors:
Marco Franzon <marco.franzon@exact-lab.it>
Tommaso Comparin <tommaso.comparin@exact-lab.it>
Jacopo Nespolo <jacopo.nespolo@exact-lab.it>


This file is part of Fractal and was originally developed by eXact lab S.r.l.
<exact-lab.it> under contract with Liberali Lab from the Friedrich Miescher
Institute for Biomedical Research and Pelkmans Lab from the University of
Zurich.
"""
import asyncio
from functools import partial
from functools import wraps
from typing import Callable


def async_wrap(func: Callable) -> Callable:
    """
    See issue #140 and https://stackoverflow.com/q/43241221/19085332

    By replacing
        .. = final_metadata.result()
    with
        .. = await async_wrap(get_app_future_result)(app_future=final_metadata)
    we avoid a (long) blocking statement.
    """

    @wraps(func)
    async def run(*args, loop=None, executor=None, **kwargs):
        if loop is None:
            loop = asyncio.get_event_loop()
        pfunc = partial(func, *args, **kwargs)
        return await loop.run_in_executor(executor, pfunc)

    return run
