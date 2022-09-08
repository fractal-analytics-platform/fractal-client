#!/bin/bash

ROOTDIR="../../.."
HERE=`pwd`

cd $ROOTDIR

rm -r data
rm -r fractal_server/migrations/versions
mkdir data
mkdir fractal_server/migrations/versions

echo -e "DEPLOYMENT_TYPE=testing\nJWT_SECRET_KEY=secret\nDATA_DIR_ROOT=/tmp/\nSQLITE_PATH=./data/fractal_server.db" > .fractal_server.env

poetry run alembic revision --autogenerate -m 'init'
poetry run alembic upgrade head

cd $HERE

echo -e "DEPLOYMENT_TYPE=testing\nJWT_SECRET_KEY=secret\nDATA_DIR_ROOT=/tmp/\nSQLITE_PATH=./data/fractal_server.db" > .fractal_server.env

rm -r data
mv $ROOTDIR/data .
