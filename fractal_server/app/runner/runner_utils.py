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
import logging
from functools import partial
from functools import wraps
from typing import Callable

import parsl
from parsl.addresses import address_by_hostname
from parsl.channels import LocalChannel
from parsl.config import Config
from parsl.dataflow.dflow import DataFlowKernelLoader
from parsl.executors import HighThroughputExecutor
from parsl.launchers import SingleNodeLauncher
from parsl.launchers import SrunLauncher
from parsl.monitoring.monitoring import MonitoringHub
from parsl.providers import LocalProvider
from parsl.providers import SlurmProvider

from ...config import settings


def add_prefix(*, workflow_id: int, executor_label: str):
    return f"{workflow_id}___{executor_label}"


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


def load_parsl_config(
    *,
    workflow_id: int,
    enable_monitoring: bool = True,
    logger: logging.Logger = None,
):

    if logger is None:
        logger = logging.getLogger("logs")

    logger.info(f"{settings.PARSL_CONFIG=}")
    logger.info(f"{settings.PARSL_DEFAULT_EXECUTOR=}")

    allowed_configs = ["local", "pelkmanslab", "fmi", "custom"]
    config = settings.PARSL_CONFIG
    if config not in allowed_configs:
        raise ValueError(f"{config=} not in {allowed_configs=}")
    if config == "custom":
        raise NotImplementedError

    if config == "local":
        # Define a single provider
        prov_local = LocalProvider(
            launcher=SingleNodeLauncher(debug=False),
            channel=LocalChannel(),
            init_blocks=1,
            min_blocks=0,
            max_blocks=4,
        )
        # Define executors
        providers = [prov_local] * 4
        labels = ["cpu-low", "cpu-mid", "cpu-high", "gpu"]
        executors = []
        for provider, label in zip(providers, labels):
            executors.append(
                HighThroughputExecutor(
                    label=add_prefix(
                        workflow_id=workflow_id, executor_label=label
                    ),
                    provider=provider,
                    address=address_by_hostname(),
                    cpu_affinity="block",
                )
            )

    elif config == "pelkmanslab":

        # Define providers
        common_args = dict(
            partition="main",
            nodes_per_block=1,
            init_blocks=1,
            min_blocks=0,
            max_blocks=100,
            parallelism=1,
            exclusive=False,
            walltime="20:00:00",
        )
        prov_cpu_low = SlurmProvider(
            launcher=SrunLauncher(debug=False),
            channel=LocalChannel(),
            cores_per_node=1,
            mem_per_node=7,
            **common_args,
        )
        prov_cpu_mid = SlurmProvider(
            launcher=SrunLauncher(debug=False),
            channel=LocalChannel(),
            cores_per_node=4,
            mem_per_node=15,
            **common_args,
        )
        prov_cpu_high = SlurmProvider(
            launcher=SrunLauncher(debug=False),
            channel=LocalChannel(),
            cores_per_node=16,
            mem_per_node=63,
            **common_args,
        )
        # Define a gpu provider
        prov_gpu = SlurmProvider(
            partition="gpu",
            launcher=SrunLauncher(debug=False),
            channel=LocalChannel(),
            nodes_per_block=1,
            init_blocks=1,
            min_blocks=0,
            max_blocks=2,
            mem_per_node=61,
            walltime="10:00:00",
        )

        # Define executors
        providers = [prov_cpu_low, prov_cpu_mid, prov_cpu_high, prov_gpu]
        labels = ["cpu-low", "cpu-mid", "cpu-high", "gpu"]
        # FIXME
        list_mem_per_worker = [7, 15, 63, 61]  # FIXME
        executors = []
        for provider, label in zip(providers, labels):
            executors.append(
                HighThroughputExecutor(
                    label=add_prefix(
                        workflow_id=workflow_id, executor_label=label
                    ),
                    provider=provider,
                    mem_per_worker=list_mem_per_worker[
                        labels.index(label)
                    ],  # FIXME
                    max_workers=100,
                    address=address_by_hostname(),
                    cpu_affinity="block",
                )
            )

    elif config == "fmi":

        # FIXME: fix partition

        # Define a cpu provider
        provider_args = dict(
            partition="main",
            launcher=SrunLauncher(debug=False),
            channel=LocalChannel(),
            nodes_per_block=1,
            init_blocks=1,
            min_blocks=0,
            max_blocks=4,
            walltime="10:00:00",
            exclusive=False,
        )
        prov_slurm_cpu = SlurmProvider(**provider_args)

        # Define executors
        htex_slurm_cpu = HighThroughputExecutor(
            label=add_prefix(
                workflow_id=workflow_id, executor_label="cpu-low"
            ),
            provider=prov_slurm_cpu,
            address=address_by_hostname(),
            cpu_affinity="block",
        )
        htex_slurm_cpu_2 = HighThroughputExecutor(
            label=add_prefix(workflow_id=workflow_id, executor_label="cpu-2"),
            provider=prov_slurm_cpu,
            address=address_by_hostname(),
            cpu_affinity="block",
        )

        executors = [htex_slurm_cpu, htex_slurm_cpu_2]

    # Extract the executor labels
    new_executor_labels = [executor.label for executor in executors]

    # Define monitoring hub and finalize configuration
    if enable_monitoring:
        monitoring = MonitoringHub(
            hub_address=address_by_hostname(),
            workflow_name="fractal",
        )
    else:
        monitoring = None

    try:
        dfk = DataFlowKernelLoader.dfk()
        old_executor_labels = [
            executor_label for executor_label in dfk.executors.keys()
        ]
        logger.info(
            f"DFK {dfk} exists, with {len(dfk.executors)} executors: "
            f"{old_executor_labels}"
        )
        logger.info(
            f"Adding {len(executors)} new executors: {new_executor_labels}"
        )

        # FIXME: what if an executor was already there?
        # (re-submitting same workflow?)

        dfk.add_executors(executors)

    # FIXME: better exception handling
    except RuntimeError:
        config = Config(
            executors=executors, monitoring=monitoring, max_idletime=20.0
        )
        logger.info(
            "DFK probably missing, "
            "proceed with parsl.clear and parsl.config.Config"
        )
        parsl.clear()
        parsl.load(config)
    dfk = DataFlowKernelLoader.dfk()
    executor_labels = [
        executor_label for executor_label in dfk.executors.keys()
    ]
    logger.info(
        f"DFK {dfk} now has {len(executor_labels)} executors: "
        f"{executor_labels}"
    )


def shutdown_executors(*, workflow_id: str, logger: logging.Logger = None):

    if logger is None:
        logger = logging.getLogger("logs")

    # Remove executors from parsl DFK
    # FIXME decorate with monitoring logs, as in:
    # https://github.com/Parsl/parsl/blob/master/parsl/dataflow/dflow.py#L1106
    dfk = DataFlowKernelLoader.dfk()
    for label, executor in dfk.executors.items():
        if label.startswith(f"{workflow_id}___"):
            executor.shutdown()
            logger.info(f"SHUTTING DOWN {label}")
    executor_labels = [
        executor_label for executor_label in dfk.executors.keys()
    ]
    logger.info(
        f"DFK {dfk} now has {len(executor_labels)} executors: "
        f"{executor_labels}"
    )


def get_unique_executor(*, workflow_id: int, task_executor: str = None):

    # Handle missing value
    if task_executor is None:
        task_executor = settings.PARSL_DEFAULT_EXECUTOR

    # Redefine task_executor, by prepending workflow_id
    new_task_executor = add_prefix(
        workflow_id=workflow_id, executor_label=task_executor
    )

    # Verify match between new_task_executor and available executors
    valid_executor_labels = DataFlowKernelLoader.dfk().executors.keys()
    if new_task_executor not in valid_executor_labels:
        raise ValueError(
            f"Executor label {new_task_executor} is not in "
            f"{valid_executor_labels=}"
        )

    return new_task_executor
