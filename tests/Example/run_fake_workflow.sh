#!/bin/bash

# 1 well, 2x2 sites
PATH_INPUT=/data/homes/fractal/mwe_fractal/tests/data/png
WFPARAMS=wf_params.json

PWD=`pwd`
PATH_OUTPUT=${PWD}/Output

CMD='poetry run python /data/homes/fractal/mwe_fractal/fractal/fractal_cmd.py'

# Clean up
rm -rf $PATH_OUTPUT
mkdir -p $PATH_OUTPUT

# Create project
$CMD project new mwe-test $PATH_OUTPUT dstest

# Add dataset to project
$CMD dataset add-resources mwe-test dstest $PATH_INPUT
$CMD dataset update-type mwe-test dstest png
$CMD dataset list mwe-test
echo

# Add tasks to project
$CMD task add create_zarr_structure png zarr none
$CMD task add yokogawa_to_zarr zarr zarr well
$CMD task list
echo

# Create workflow
$CMD workflow new mwe-test wftest create_zarr_structure
$CMD workflow add-task mwe-test wftest yokogawa_to_zarr
$CMD workflow list mwe-test
echo

# Execute workflow
$CMD workflow apply mwe-test wftest dstest dstest $PATH_INPUT $PATH_OUTPUT $WFPARAMS
echo
