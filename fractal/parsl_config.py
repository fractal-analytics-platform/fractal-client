from parsl.addresses import address_by_hostname
from parsl.channels import LocalChannel
from parsl.executors import HighThroughputExecutor
from parsl.launchers import SrunLauncher
from parsl.monitoring.monitoring import MonitoringHub
from parsl.providers import SlurmProvider


def define_MonitoringHub():
    return MonitoringHub(
        hub_address=address_by_hostname(),
        hub_port=55055,
        monitoring_debug=True,
        resource_monitoring_interval=5,
    )


def define_SlurmProvider(
    nodes_per_block=1,
    cores_per_node=16,
    mem_per_node_GB=10,
    partition="main",
    worker_init="",
):
    slurm = SlurmProvider(
        channel=LocalChannel(),
        nodes_per_block=nodes_per_block,
        cores_per_node=cores_per_node,
        mem_per_node=mem_per_node_GB,  # specified in GB
        partition=partition,
        worker_init=worker_init,
        launcher=SrunLauncher(debug=True),
        parallelism=0,
        walltime="10:00:00",
        cmd_timeout=60,
    )
    return slurm


def define_HighThroughputExecutor(provider=None, max_workers=40):
    htex = HighThroughputExecutor(
        label="htex",
        address=address_by_hostname(),
        worker_debug=True,
        max_workers=max_workers,
        provider=provider,
    )
    return htex
