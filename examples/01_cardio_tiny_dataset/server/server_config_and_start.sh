#!/bin/bash

# Start with an empty db (created from within the fractal-server repo)
rm -r data
mkdir data
cp -v empty_db/fractal_server.db data/

# Set environment variables
echo -e "DEPLOYMENT_TYPE=testing\nJWT_SECRET_KEY=secret\nDATA_DIR_ROOT=/tmp/\nSQLITE_PATH=./data/fractal_server.db" > .fractal_server.env

# Start the server (note that this may be renamed later, see https://github.com/fractal-analytics-platform/fractal-server/issues/21)
server
