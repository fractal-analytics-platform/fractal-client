PATH_INPUT=/data/active/fractal-temp/3D/PelkmansLab/CardiacMultiplexing/Cycle1_testSubset/
#PATH_INPUT=/data/active/fractal-temp/3D/PelkmansLab/CardiacMultiplexing/Cycle1_subset/

AUXDIR=/data/homes/fractal/mwe_fractal_tommaso/fractal/Temporary_folder_for_regression_test
PATH_OUTPUT=${AUXDIR}/result/
CMD=fractal_cmd.py
CMD=fractal_cmd_parsl.py

#echo 'Re-install poetry'
#poetry install

echo 'Clean up old dstest folder'
rm -rf $AUXDIR
mkdir -p $AUXDIR

echo 'Create project'
poetry run python $CMD project new mwe-test $AUXDIR dstest
echo

echo 'Add dataset to project'
poetry run python $CMD dataset add-resources mwe-test dstest $PATH_INPUT
poetry run python $CMD dataset list mwe-test
echo

echo 'Update dataset type'
poetry run python $CMD dataset update-type mwe-test dstest png
poetry run python $CMD dataset list mwe-test
echo

echo 'Add create_zarr_structure task'
poetry run python $CMD task add create_zarr_structure png zarr
poetry run python $CMD task list
echo

echo 'Add yokogawa_to_zarr task'
poetry run python $CMD task add yokogawa_to_zarr zarr zarr
poetry run python $CMD task list
echo


echo 'Add maximum_intensity_projection'
poetry run python $CMD task add maximum_intensity_projection zarr zarr
poetry run python $CMD task list
echo

echo 'Create workflow'
poetry run python $CMD workflow new mwe-test wftest create_zarr_structure
poetry run python $CMD workflow list mwe-test
echo

echo 'Add yokogawa_to_zarr task'
poetry run python $CMD workflow add-task mwe-test wftest yokogawa_to_zarr
poetry run python $CMD workflow list mwe-test
echo

echo 'Add maximum_intensity_projection'
poetry run python $CMD workflow add-task mwe-test wftest maximum_intensity_projection
poetry run python $CMD workflow list mwe-test
echo

echo 'Execute workflow'
poetry run python $CMD workflow apply mwe-test wftest dstest dstest $PATH_INPUT $PATH_OUTPUT
echo
