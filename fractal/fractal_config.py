# Command to be included at the beginning of SLURM submission scripts
worker_init = "source /opt/easybuild/software/Anaconda3/2019.07/"
worker_init += "etc/profile.d/conda.sh\n"
worker_init += "conda activate fractal"

max_workers = 32

# SLURM options to be passed to SlurmProvider in parsl
nodes_per_block = 1
cores_per_node = 16
mem_per_node_GB = 64
partition = "main"
