"""
Copyright 2022 (C) Friedrich Miescher Institute for Biomedical Research and
University of Zurich

Original authors:
Tommaso Comparin <tommaso.comparin@exact-lab.it>
Marco Franzon <marco.franzon@exact-lab.it>

This file is part of Fractal and was originally developed by eXact lab S.r.l.
<exact-lab.it> under contract with Liberali Lab from the Friedrich Miescher
Institute for Biomedical Research and Pelkmans Lab from the University of
Zurich.
"""
from parsl.addresses import address_by_hostname
from parsl.channels import LocalChannel
from parsl.executors import HighThroughputExecutor
from parsl.launchers import SrunLauncher
from parsl.monitoring.monitoring import MonitoringHub
from parsl.providers import SlurmProvider


def define_MonitoringHub(workflow_name=None):
    # FIXME: scan ports and find one that is available
    kwargs = dict(
        hub_address=address_by_hostname(),
        # monitoring_debug=True,
        resource_monitoring_interval=30,
    )
    if workflow_name is not None:
        kwargs["workflow_name"] = workflow_name
    return MonitoringHub(**kwargs)


def define_SlurmProvider(
    nodes_per_block=1,
    cores_per_node=16,
    mem_per_node_GB=10,
    partition="main",
    worker_init="",
    max_blocks=1,
    exclusive=False,
):

    slurm = SlurmProvider(
        channel=LocalChannel(),
        nodes_per_block=nodes_per_block,
        cores_per_node=cores_per_node,
        mem_per_node=mem_per_node_GB,  # specified in GB
        partition=partition,
        worker_init=worker_init,
        launcher=SrunLauncher(),  # debug=True),
        walltime="23:00:00",
        cmd_timeout=60,
        min_blocks=1,
        max_blocks=max_blocks,
        parallelism=1,
        exclusive=exclusive,
    )

    return slurm


def define_HighThroughputExecutor(provider=None, max_workers=40, label="htex"):

    htex = HighThroughputExecutor(
        label=label,
        address=address_by_hostname(),
        # worker_debug=True,
        max_workers=max_workers,
        cores_per_worker=16,
        provider=provider,
        cpu_affinity="block",
    )
    return htex
