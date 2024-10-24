#!/bin/sh

export DB_ENGINE=postgres-psycopg
export POSTGRES_HOST=localhost
export POSTGRES_DB=pytest-fractal-client
export POSTGRES_USER=postgres
export POSTGRES_PASSWORD=postgres
export FRACTAL_RUNNER_BACKEND=local
export JWT_SECRET_KEY=secret_key
export FRACTAL_TASKS_DIR=FRACTAL_TASK_DIR
export FRACTAL_RUNNER_WORKING_BASE_DIR=FRACTAL_RUNNER_WORKING_BASE_DIR
export FRACTAL_LOGGING_LEVEL=0

poetry run fractalctl set-db
poetry run fractalctl start --port 8765 > /dev/null 2>&1 &
