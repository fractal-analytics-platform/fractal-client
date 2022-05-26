# Parameters of parsl.executors.HighThroughputExecutor
max_workers = 32

# Parameters of parsl.providers.SlurmProvider
# Note that worker_init is a command which is included at the beginning of
# each SLURM submission scripts
nodes_per_block = 1
cores_per_node = 16
mem_per_node_GB = 64
partition = "main"
worker_init = "source /opt/easybuild/software/Anaconda3/2019.07/"
worker_init += "etc/profile.d/conda.sh\n"
worker_init += "conda activate fractal"
