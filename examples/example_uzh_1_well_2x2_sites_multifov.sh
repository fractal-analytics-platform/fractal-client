#!/bin/bash

#2x2, single well
PATH_INPUT=/data/active/fractal/3D/PelkmansLab/CardiacMultiplexing/Cycle1_testSubset
WFPARAMS=wf_params_uzh_1_well_2x2_sites.json

MWE_DIR=/data/active/fractal/tests
PATH_OUTPUT=${MWE_DIR}/Temporary_data_UZH_1_well_2x2_sites_multifov

CMD='poetry run python ../fractal/fractal_cmd.py'

#echo 'Re-install poetry'
#poetry install

echo 'Clean up'
rm -rf $PATH_OUTPUT
mkdir -p $PATH_OUTPUT

echo 'Create project'
$CMD project new mwe-test $PATH_OUTPUT dstest
echo

echo 'Add dataset to project'
$CMD dataset add-resources mwe-test dstest $PATH_INPUT
$CMD dataset list mwe-test
echo

echo 'Update dataset type'
$CMD dataset update-type mwe-test dstest png
$CMD dataset list mwe-test
echo

echo 'Add create_zarr_structure_multifov task'
$CMD task add create_zarr_structure_multifov png zarr none
$CMD task list
echo

echo 'Add yokogawa_to_zarr_multifov task'
$CMD task add yokogawa_to_zarr_multifov zarr zarr well
$CMD task list
echo

echo 'Create workflow'
$CMD workflow new mwe-test wftest create_zarr_structure_multifov
$CMD workflow list mwe-test
echo

echo 'Add yokogawa_to_zarr_multifov task'
$CMD workflow add-task mwe-test wftest yokogawa_to_zarr_multifov
$CMD workflow list mwe-test
echo

echo 'Execute workflow'
$CMD workflow apply mwe-test wftest dstest dstest $PATH_INPUT $PATH_OUTPUT $WFPARAMS
echo
