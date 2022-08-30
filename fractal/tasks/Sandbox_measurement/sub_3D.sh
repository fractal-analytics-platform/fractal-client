#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --partition=main
#SBATCH --time=10:00:00

date
echo

poetry run python call_measurement_3D.py

echo
date
