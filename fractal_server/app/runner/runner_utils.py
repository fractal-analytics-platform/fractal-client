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
import os
from functools import partial
from functools import wraps
from typing import Callable

from parsl.addresses import address_by_hostname
from parsl.channels import LocalChannel
from parsl.config import Config
from parsl.dataflow.dflow import DataFlowKernel
from parsl.executors import HighThroughputExecutor
from parsl.launchers import SingleNodeLauncher
from parsl.launchers import SrunLauncher
from parsl.monitoring.monitoring import MonitoringHub
from parsl.providers import LocalProvider
from parsl.providers import SlurmProvider

from ...config_runner import settings


def add_prefix(*, workflow_id: int, executor_label: str):
    return f"{workflow_id}__{executor_label}"


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


class LocalChannel_fractal(LocalChannel):
    def __init__(self, *args, username: str = None, **kwargs):
        self.username: str = username
        super().__init__(*args, **kwargs)

    def makedirs(self, path, mode=0o777, exist_ok=False):
        """
        Changes with respect to LocalChannel.makedirs:
          * Set mode=0o777 default
          * Set umask(0), to mask user limitations, and reset it after mkdir.
        """
        old_umask = os.umask(0)
        os.makedirs(path, mode, exist_ok)
        os.umask(old_umask)
        return f"Create {path} dir"

    def execute_wait(self, cmd, *args, **kwargs):
        if self.username is None:
            return super().execute_wait(cmd, *args, **kwargs)
        else:
            if cmd.startswith(("/bin/bash -c '", "ps", "kill")):
                msg = (
                    f'We cannot add "sudo - {self.username} -c" in front of '
                    f"LocalProvider commands like {cmd=}"
                )
                raise NotImplementedError(msg)
            elif cmd.startswith(("sbatch", "scancel")):
                if cmd == "scancel ":
                    new_cmd = (
                        'echo "execute_wait was called with dummy "'
                        f'"command "{cmd}". We are replacing it with '
                        'this echo command."'
                    )
                else:
                    new_cmd = f'sudo su - {self.username} -c "{cmd}"'
            elif cmd.startswith("squeue"):
                new_cmd = cmd
            else:
                msg = "We cannot add sudo in front of " f"commands like {cmd=}"
                raise NotImplementedError(msg)
            return super().execute_wait(new_cmd, *args, **kwargs)


def generate_parsl_config(
    *,
    workflow_id: int,
    workflow_name: str,
    enable_monitoring: bool = True,
    username: str = None,
    worker_init: str = None,
) -> Config:

    allowed_configs = ["local", "pelkmanslab", "fmi", "custom"]
    config = settings.RUNNER_CONFIG
    if config not in allowed_configs:
        raise ValueError(f"{config=} not in {allowed_configs=}")
    if config == "custom":
        raise NotImplementedError

    # Prepare log folders for channels and executors
    RUNNER_LOG_DIR = settings.RUNNER_LOG_DIR
    if not os.path.isdir(RUNNER_LOG_DIR):
        old_umask = os.umask(0)
        os.mkdir(RUNNER_LOG_DIR)
        os.umask(old_umask)
    workflow_log_dir = f"{RUNNER_LOG_DIR}/workflow_{workflow_id:06d}"
    workflow_log_dir = os.path.abspath(workflow_log_dir)
    if not os.path.isdir(workflow_log_dir):
        old_umask = os.umask(0)
        os.mkdir(workflow_log_dir, mode=0o777)
        os.umask(old_umask)
    script_dir = f"{workflow_log_dir}/scripts"
    script_dir = os.path.abspath(script_dir)
    if not os.path.isdir(script_dir):
        old_umask = os.umask(0)
        os.mkdir(script_dir)
        os.umask(old_umask)
    worker_logdir_root = f"{workflow_log_dir}/executors"
    if not os.path.isdir(worker_logdir_root):
        old_umask = os.umask(0)
        os.mkdir(worker_logdir_root, 0o777)
        os.umask(old_umask)
    channel_args = dict(script_dir=script_dir)

    # Set additional parameters for channel/provider
    if username is not None:
        channel_args["username"] = username
    if worker_init is None:
        worker_init = ""

    if config == "local":
        # Define a single provider
        prov_local = LocalProvider(
            move_files=True,
            launcher=SingleNodeLauncher(debug=False),
            channel=LocalChannel_fractal(**channel_args),
            init_blocks=1,
            min_blocks=0,
            max_blocks=4,
            worker_init=worker_init,
        )
        # Define executors
        labels = ["cpu-low", "cpu-mid", "cpu-high", "gpu"]
        executors = []
        for label in labels:
            htex = HighThroughputExecutor(
                label=add_prefix(
                    workflow_id=workflow_id, executor_label=label
                ),
                provider=prov_local,
                address=address_by_hostname(),
                cpu_affinity="block",
                worker_logdir_root=worker_logdir_root,
            )
            executors.append(htex)

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
            move_files=True,
            worker_init=worker_init,
        )
        prov_cpu_low = SlurmProvider(
            launcher=SrunLauncher(debug=False),
            channel=LocalChannel_fractal(**channel_args),
            cores_per_node=1,
            mem_per_node=7,
            **common_args,
        )
        prov_cpu_mid = SlurmProvider(
            launcher=SrunLauncher(debug=False),
            channel=LocalChannel_fractal(**channel_args),
            cores_per_node=4,
            mem_per_node=15,
            **common_args,
        )
        prov_cpu_high = SlurmProvider(
            launcher=SrunLauncher(debug=False),
            channel=LocalChannel_fractal(**channel_args),
            cores_per_node=16,
            mem_per_node=61,
            **common_args,
        )
        # Define a gpu provider
        prov_gpu = SlurmProvider(
            partition="gpu",
            launcher=SrunLauncher(debug=False),
            channel=LocalChannel_fractal(**channel_args),
            nodes_per_block=1,
            init_blocks=1,
            min_blocks=0,
            max_blocks=2,
            mem_per_node=61,
            walltime="20:00:00",
            move_files=True,
            worker_init=worker_init,
        )

        # Define executors
        providers = [prov_cpu_low, prov_cpu_mid, prov_cpu_high, prov_gpu]
        labels = ["cpu-low", "cpu-mid", "cpu-high", "gpu"]
        executors = []
        for provider, label in zip(providers, labels):
            htex = HighThroughputExecutor(
                label=add_prefix(
                    workflow_id=workflow_id, executor_label=label
                ),
                provider=provider,
                mem_per_worker=provider.mem_per_node,
                max_workers=100,
                address=address_by_hostname(),
                cpu_affinity="block",
                worker_logdir_root=worker_logdir_root,
            )
            executors.append(htex)

    elif config == "fmi":
        common_args = dict(
            partition="cpu_long",
            nodes_per_block=1,
            init_blocks=1,
            min_blocks=0,
            max_blocks=100,
            parallelism=1,
            exclusive=False,
            walltime="20:00:00",
            move_files=True,
            worker_init=worker_init,
        )
        prov_cpu_low = SlurmProvider(
            launcher=SrunLauncher(debug=False),
            channel=LocalChannel_fractal(**channel_args),
            cores_per_node=1,
            mem_per_node=7,
            **common_args,
        )
        prov_cpu_mid = SlurmProvider(
            launcher=SrunLauncher(debug=False),
            channel=LocalChannel_fractal(**channel_args),
            cores_per_node=4,
            mem_per_node=15,
            **common_args,
        )
        prov_cpu_high = SlurmProvider(
            launcher=SrunLauncher(debug=False),
            channel=LocalChannel_fractal(**channel_args),
            cores_per_node=20,
            mem_per_node=61,
            **common_args,
        )

        # Define executors
        providers = [prov_cpu_low, prov_cpu_mid, prov_cpu_high]
        labels = ["cpu-low", "cpu-mid", "cpu-high"]

        executors = []
        for provider, label in zip(providers, labels):
            executors.append(
                HighThroughputExecutor(
                    label=add_prefix(
                        workflow_id=workflow_id, executor_label=label
                    ),
                    provider=provider,
                    mem_per_worker=provider.mem_per_node,
                    max_workers=100,
                    address=address_by_hostname(),
                    cpu_affinity="block",
                    worker_logdir_root=worker_logdir_root,
                )
            )

    # Define monitoring hub and finalize configuration
    if enable_monitoring:
        monitoring = MonitoringHub(
            hub_address=address_by_hostname(),
            workflow_name=workflow_name,
        )
    else:
        monitoring = None

    config = Config(
        executors=executors, monitoring=monitoring, max_idletime=20.0
    )
    return config


def load_parsl_config(
    *,
    workflow_id: int,
    workflow_name: str = "default workflow name",
    config: Config = None,
    enable_monitoring: bool = True,
    username: str = None,
    worker_init: str = None,
):

    logger = logging.getLogger(f"WF{workflow_id}")
    logger.info(f"{settings.RUNNER_CONFIG=}")
    logger.info(f"{settings.RUNNER_DEFAULT_EXECUTOR=}")

    if not config:
        config = generate_parsl_config(
            enable_monitoring=enable_monitoring,
            workflow_id=workflow_id,
            workflow_name=workflow_name,
            username=username,
            worker_init=worker_init,
        )

    """
    try:
        dfk = DataFlowKernelLoader.dfk()
        old_executor_labels = [
            executor_label for executor_label in dfk.executors.keys()
        ]
        logger.info(
            f"DFK {dfk} exists, with {len(dfk.executors)} executors: "
            f"{old_executor_labels}"
        )

        # FIXME: what if an executor was already there?
        # (re-submitting same workflow?)
        dfk.add_executors(config.executors)

    # FIXME: better exception handling
    except RuntimeError:
        logger.info(
            "DFK probably missing, "
            "proceed with parsl.clear and parsl.config.Config"
        )
        parsl.clear()
        DataFlowKernelLoader.load(config)
        dfk = DataFlowKernelLoader.dfk()
    """

    dfk = DataFlowKernel(config=config)

    executor_labels = [
        executor_label for executor_label in dfk.executors.keys()
    ]
    logger.info(
        f"DFK {dfk} now has {len(executor_labels)} executors: "
        f"{executor_labels}"
    )

    return dfk


def get_unique_executor(
    *, workflow_id: int, task_executor: str = None, data_flow_kernel=None
):

    # Handle missing value
    if task_executor is None:
        task_executor = settings.RUNNER_DEFAULT_EXECUTOR

    # Redefine task_executor, by prepending workflow_id
    new_task_executor = add_prefix(
        workflow_id=workflow_id, executor_label=task_executor
    )

    # Verify match between new_task_executor and available executors
    valid_executor_labels = data_flow_kernel.executors.keys()
    if new_task_executor not in valid_executor_labels:
        raise ValueError(
            f"Executor label {new_task_executor} is not in "
            f"{valid_executor_labels=}"
        )

    return new_task_executor
