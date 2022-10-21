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
import logging
import os
import warnings
from pathlib import Path

from parsl.addresses import address_by_hostname
from parsl.channels import LocalChannel
from parsl.config import Config as ParslConfig
from parsl.dataflow.dflow import DataFlowKernel
from parsl.executors import HighThroughputExecutor
from parsl.launchers import SingleNodeLauncher
from parsl.launchers import SrunLauncher
from parsl.monitoring.monitoring import MonitoringHub
from parsl.providers import LocalProvider
from parsl.providers import SlurmProvider

from ....config_runner import settings


def get_executor_label(*, workflow_id: int, executor_label: str):
    return f"{workflow_id}__{executor_label}"


class FractalLocalChannel(LocalChannel):
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
                    f'We do not add "sudo - {self.username} -c" in front of '
                    f"LocalProvider commands like {cmd=}."
                )
                warnings.warn(msg)
                new_cmd = cmd
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
    workflow_dir: Path,
    enable_monitoring: bool = True,
    username: str = None,
    worker_init: str = None,
) -> ParslConfig:

    if enable_monitoring:
        msg = (
            "parsl monitoring is currently disabled, due to"
            "https://github.com/fractal-analytics-platform/"
            "fractal-server/issues/148"
        )
        raise NotImplementedError(msg)

    allowed_configs = ["minimal", "local", "pelkmanslab", "fmi", "custom"]
    config = settings.RUNNER_CONFIG
    if config not in allowed_configs:
        raise ValueError(f"{config=} not in {allowed_configs=}")
    if config == "custom":
        raise NotImplementedError

    # Prepare log folders for channels and executors
    script_dir = workflow_dir / "scripts"
    worker_logdir_root = workflow_dir / "executors"

    script_dir.mkdir(exist_ok=True, parents=True)
    worker_logdir_root.mkdir(exist_ok=True, parents=True)

    channel_args = dict(script_dir=script_dir)

    # Set additional parameters for channel/provider
    if username is not None:
        channel_args["username"] = username
    if worker_init is None:
        worker_init = ""

    if config == "minimal":
        # Define a single provider
        prov_local = LocalProvider(
            move_files=True,
            launcher=SingleNodeLauncher(debug=False),
            channel=FractalLocalChannel(**channel_args),
            init_blocks=1,
            min_blocks=0,
            max_blocks=4,
            worker_init=worker_init,
        )
        # Define executor
        executors = [
            HighThroughputExecutor(
                label=get_executor_label(
                    workflow_id=workflow_id, executor_label="cpu-low"
                ),
                provider=prov_local,
                address=address_by_hostname(),
                cpu_affinity="block",
                worker_logdir_root=worker_logdir_root,
            )
        ]

    elif config == "local":
        # Define a single provider
        prov_local = LocalProvider(
            move_files=True,
            launcher=SingleNodeLauncher(debug=False),
            channel=FractalLocalChannel(**channel_args),
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
                label=get_executor_label(
                    workflow_id=workflow_id, executor_label=label
                ),
                provider=prov_local,
                address=address_by_hostname(),
                cpu_affinity="block",
                worker_logdir_root=worker_logdir_root.as_posix(),
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
            channel=FractalLocalChannel(**channel_args),
            cores_per_node=1,
            mem_per_node=7,
            **common_args,
        )
        prov_cpu_mid = SlurmProvider(
            launcher=SrunLauncher(debug=False),
            channel=FractalLocalChannel(**channel_args),
            cores_per_node=4,
            mem_per_node=15,
            **common_args,
        )
        prov_cpu_high = SlurmProvider(
            launcher=SrunLauncher(debug=False),
            channel=FractalLocalChannel(**channel_args),
            cores_per_node=16,
            mem_per_node=61,
            **common_args,
        )
        # Define a gpu provider
        prov_gpu = SlurmProvider(
            partition="gpu",
            launcher=SrunLauncher(debug=False),
            channel=FractalLocalChannel(**channel_args),
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
                label=get_executor_label(
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
            channel=FractalLocalChannel(**channel_args),
            cores_per_node=1,
            mem_per_node=7,
            **common_args,
        )
        prov_cpu_mid = SlurmProvider(
            launcher=SrunLauncher(debug=False),
            channel=FractalLocalChannel(**channel_args),
            cores_per_node=4,
            mem_per_node=15,
            **common_args,
        )
        prov_cpu_high = SlurmProvider(
            launcher=SrunLauncher(debug=False),
            channel=FractalLocalChannel(**channel_args),
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
                    label=get_executor_label(
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

    config = ParslConfig(
        executors=executors, monitoring=monitoring, max_idletime=20.0
    )
    return config


def load_parsl_config(
    *,
    workflow_id: int,
    logger: logging.Logger,
    workflow_name: str,
    workflow_dir: Path,
    parsl_config: ParslConfig = None,
    enable_monitoring: bool = None,
    username: str = None,
    worker_init: str = None,
):

    if enable_monitoring is None:
        enable_monitoring = settings.RUNNER_MONITORING

    if not parsl_config:
        parsl_config = generate_parsl_config(
            enable_monitoring=enable_monitoring,
            workflow_id=workflow_id,
            workflow_name=workflow_name,
            workflow_dir=workflow_dir,
            username=username,
            worker_init=worker_init,
        )

    dfk = DataFlowKernel(config=parsl_config)

    executor_labels = [
        executor_label for executor_label in dfk.executors.keys()
    ]
    logger.info(f"New DFK {dfk}, with executors {executor_labels}")

    return dfk
