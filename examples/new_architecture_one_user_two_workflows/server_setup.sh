#!/bin/bash

ROOTDIR="../.."
HERE=`pwd`

cd $ROOTDIR

rm -r data
rm -r fractal/server/migrations/versions
mkdir data
mkdir fractal/server/migrations/versions

echo -e "DEPLOYMENT_TYPE=testing\nJWT_SECRET_KEY=secret\nDATA_DIR_ROOT=/tmp/\nSQLITE_PATH=./data/fractal_server.db" > .fractal_server.env
echo -e "FRACTAL_USER=test@me.com\nFRACTAL_PASSWORD=test" > .fractal.env

poetry run alembic revision --autogenerate -m 'init'
poetry run alembic upgrade head

cd $HERE

echo -e "DEPLOYMENT_TYPE=testing\nJWT_SECRET_KEY=secret\nDATA_DIR_ROOT=/tmp/\nSQLITE_PATH=./data/fractal_server.db" > .fractal_server.env
echo -e "FRACTAL_USER=test@me.com\nFRACTAL_PASSWORD=test" > .fractal.env

rm -r data
mv $ROOTDIR/data .
