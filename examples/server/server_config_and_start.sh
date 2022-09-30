#!/bin/bash

# Create an empty db
rm fractal_server/migrations/versions/*.py -v
rm -r data
mkdir data
alembic revision --autogenerate
alembic upgrade head

# Remove old stuff
# rm -r runinfo
# rm fractal.log
rm cmd_parsl.slurm.*.*.sh
rm -v parsl_executors.log

# Set environment variables
echo -e "DEPLOYMENT_TYPE=testing\nJWT_SECRET_KEY=secret\nDATA_DIR_ROOT=/tmp/\nSQLITE_PATH=./data/fractal_server.db" > .fractal_server.env

# Start the server
fractal-server
