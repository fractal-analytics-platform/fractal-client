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
import sys
import warnings
from contextlib import contextmanager
from pathlib import Path
from typing import Any
from typing import Dict
from typing import Optional

import parsl.executors.high_throughput
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

from ....config import get_settings
from ....syringe import Inject


class FractalLocalChannel(LocalChannel):
    def __init__(self, *args, username: str = None, **kwargs):
        self.username: Optional[str] = username
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
                    f"The attribute username={self.username} is explicitly set"
                    " but will be ignored, since we are currently using parsl "
                    "LocalProvider."
                )
                warnings.warn(msg)
                new_cmd = cmd
            elif cmd.startswith(("sbatch", "scancel")):
                new_cmd = f"sudo --non-interactive -u {self.username} {cmd}"
            elif cmd.startswith("squeue"):
                new_cmd = cmd
            else:
                msg = "We cannot add sudo in front of commands like {cmd=}"
                raise NotImplementedError(msg)
            return super().execute_wait(new_cmd, *args, **kwargs)


class FractalHighThroughputExecutor(HighThroughputExecutor):
    """
    This subclass is used to inject an executable process_worker_pool.py
    command (made by an absolute path to a python interpreter and an absolute
    path to process_worker_pool.py). This is necessary since
    process_worker_pool.py may have to be run by a different user, e.g. on a
    SLURM cluster.
    """

    def __init__(self, *args, logger_name=None, **kwargs):

        super().__init__(*args, **kwargs)

        # Assemble command for process_worker_pool.py script
        python_bin = sys.executable
        pwp_path = str(
            Path(parsl.executors.high_throughput.__file__).parent
            / "process_worker_pool.py"
        )
        pwp_command = f"{python_bin} {pwp_path}"

        # Write logs
        logger = logging.getLogger(logger_name)
        logger.info(
            "Replacing process_worker_pool.py with custom command, so that "
            "other SLURM users can execute it."
        )
        logger.warning(
            "Replacing process_worker_pool.py with custom command may "
            f"break if {python_bin} is not accessible from SLURM "
            "computing nodes."
        )

        # Replace default script with custom one
        self.launch_cmd = self.launch_cmd.replace(
            "process_worker_pool.py", pwp_command
        )


def generate_parsl_config(
    *,
    workflow_id: int,
    workflow_name: str,
    workflow_dir: Path,
    enable_monitoring: bool = True,
    username: str = None,
    worker_init: str = None,
    logger_name: str,
) -> ParslConfig:

    if enable_monitoring:
        msg = (
            "parsl monitoring is currently disabled, due to"
            "https://github.com/fractal-analytics-platform/"
            "fractal-server/issues/148"
        )
        raise NotImplementedError(msg)

    allowed_configs = ["minimal", "local", "pelkmanslab", "fmi", "custom"]
    settings = Inject(get_settings)
    config = settings.RUNNER_CONFIG
    if config not in allowed_configs:
        raise ValueError(f"{config=} not in {allowed_configs=}")
    if config == "custom":
        raise NotImplementedError

    # Prepare log folders for channels and executors
    script_dir = workflow_dir / "scripts"
    worker_logdir_root = workflow_dir / "executors"

    # Create folders with fully open permissions
    old_umask = os.umask(0)
    script_dir.mkdir(mode=0o777, exist_ok=True, parents=True)
    worker_logdir_root.mkdir(mode=0o777, exist_ok=True, parents=True)
    os.umask(old_umask)

    channel_args: Dict[str, Any] = dict(script_dir=script_dir)

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
                label="cpu-low",
                provider=prov_local,
                address=address_by_hostname(),
                cpu_affinity="block",
                worker_logdir_root=worker_logdir_root.as_posix(),
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
                label=label,
                provider=prov_local,
                address=address_by_hostname(),
                cpu_affinity="block",
                worker_logdir_root=worker_logdir_root.as_posix(),
            )
            executors.append(htex)

    elif config == "pelkmanslab":
        # Define providers
        common_args: Dict[str, Any] = dict(
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
            htex = FractalHighThroughputExecutor(
                label=label,
                provider=provider,
                mem_per_worker=provider.mem_per_node,
                max_workers=100,
                address=address_by_hostname(),
                cpu_affinity="block",
                worker_logdir_root=worker_logdir_root.as_posix(),
                logger_name=logger_name,
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
                FractalHighThroughputExecutor(
                    label=label,
                    provider=provider,
                    mem_per_worker=provider.mem_per_node,
                    max_workers=100,
                    address=address_by_hostname(),
                    cpu_affinity="block",
                    worker_logdir_root=worker_logdir_root.as_posix(),
                    logger_name=logger_name,
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

    parsl_config = ParslConfig(
        executors=executors, monitoring=monitoring, max_idletime=20.0
    )
    return parsl_config


@contextmanager
def load_parsl_config(
    *,
    workflow_id: int,
    logger_name: str,
    workflow_name: str,
    workflow_dir: Path,
    parsl_config: ParslConfig = None,
    enable_monitoring: bool = None,
    username: str = None,
    worker_init: str = None,
):
    settings = Inject(get_settings)
    logger = logging.getLogger(logger_name)

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
            logger_name=logger_name,
        )

    dfk = DataFlowKernel(config=parsl_config)

    executor_labels = [
        executor_label for executor_label in dfk.executors.keys()
    ]
    logger.info(f"New DFK {dfk}, with executors {executor_labels}")
    try:
        yield dfk
    finally:
        dfk.cleanup()
