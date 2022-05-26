#2x2, single well
PATH_INPUT=/data/active/fractal-temp/3D/PelkmansLab/CardiacMultiplexing/Cycle1_testSubset
WFPARAMS=wf_params_uzh_1_well_2x2_sites.json

MWE_DIR=/data/homes/fractal/mwe_fractal_tommaso
PATH_OUTPUT=${MWE_DIR}/examples/Temporary_data_UZH_1_well_2x2_sites

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

echo 'Add create_zarr_structure task'
$CMD task add create_zarr_structure png zarr none
$CMD task list
echo

echo 'Add yokogawa_to_zarr task'
$CMD task add yokogawa_to_zarr zarr zarr well
$CMD task list
echo

echo 'Add replicate_zarr_structure_mip'
$CMD task add replicate_zarr_structure_mip zarr zarr plate
$CMD task list
echo

echo 'Add maximum_intensity_projection'
$CMD task add maximum_intensity_projection zarr zarr well
$CMD task list
echo

echo 'Create workflow'
$CMD workflow new mwe-test wftest create_zarr_structure
$CMD workflow list mwe-test
echo

echo 'Add yokogawa_to_zarr task'
$CMD workflow add-task mwe-test wftest yokogawa_to_zarr
$CMD workflow list mwe-test
echo

echo 'Add replicate_zarr_structure_mip'
$CMD workflow add-task mwe-test wftest replicate_zarr_structure_mip
$CMD workflow list mwe-test
echo

echo 'Add maximum_intensity_projection'
$CMD workflow add-task mwe-test wftest maximum_intensity_projection
$CMD workflow list mwe-test
echo

echo 'Execute workflow'
$CMD workflow apply mwe-test wftest dstest dstest $PATH_INPUT $PATH_OUTPUT $WFPARAMS
echo
