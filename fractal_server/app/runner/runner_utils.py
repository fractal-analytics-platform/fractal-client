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

from parsl.addresses import address_by_hostname
from parsl.channels import LocalChannel
from parsl.executors import HighThroughputExecutor
from parsl.launchers import SingleNodeLauncher
from parsl.launchers import SrunLauncher
from parsl.providers import LocalProvider
from parsl.providers import SlurmProvider

from ...config import settings


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


def load_parsl_executors_config(*, config: str):

    allowed_configs = ["local", "pelkmanslab"]
    if config not in allowed_configs:
        raise ValueError(f"{config=} not in {allowed_configs=}")

    if config == "local":

        # Define a single provider
        prov = LocalProvider(
            launcher=SingleNodeLauncher(debug=False),
            channel=LocalChannel(),
            init_blocks=1,
            min_blocks=1,
            max_blocks=4,
        )

        # Define two identical (apart from the label) executors
        htex = HighThroughputExecutor(
            label="cpu",
            provider=prov,
            address=address_by_hostname(),
        )
        htex_2 = HighThroughputExecutor(
            label="cpu-2",
            provider=prov,
            address=address_by_hostname(),
        )
        executors = [htex, htex_2]

    elif config == "pelkmanslab":
        pass

        # Define a single provider
        provider_args = dict(
            partition=settings.SLURM_PARTITION_CPU,
            launcher=SrunLauncher(debug=False),
            channel=LocalChannel(),
            nodes_per_block=1,
            init_blocks=1,
            min_blocks=1,
            max_blocks=4,
            walltime="10:00:00",
        )
        prov_slurm_cpu = SlurmProvider(**provider_args)

        # Define two identical (apart from the label) executors
        htex_slurm_cpu = HighThroughputExecutor(
            label="cpu",
            provider=prov_slurm_cpu,
            address=address_by_hostname(),
            cpu_affinity="block",
        )
        htex_slurm_cpu_2 = HighThroughputExecutor(
            label="cpu-2",
            provider=prov_slurm_cpu,
            address=address_by_hostname(),
            cpu_affinity="block",
        )

        executors = [htex_slurm_cpu, htex_slurm_cpu_2]

    return executors
