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

# Parameters of parsl.executors.HighThroughputExecutor
max_workers = 1  # This is the maximum number of workers per block

# Parameters of parsl.providers.SlurmProvider
# Note that worker_init is a command which is included at the beginning of
# each SLURM submission scripts
nodes_per_block = 1  # This implies that a block corresponds to a node
max_blocks = 15  # Maximum number of blocks (=nodes) that parsl can use
exclusive = False
cores_per_node = 16
mem_per_node_GB = 60
partition = "main"
worker_init = "source /opt/easybuild/software/Anaconda3/2019.07/"
worker_init += "etc/profile.d/conda.sh\n"
worker_init += "conda activate fractal"
