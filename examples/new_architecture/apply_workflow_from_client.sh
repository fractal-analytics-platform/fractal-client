rm -r tmp-proj
mkdir tmp-proj

# Create project
poetry run client project new test-prj tmp-proj


TESTDATA=/data/homes/fractal/mwe_fractal/tests/data


PROJECT_ID=1
DATASET_ID=1
DATASET_ID_OUT=2
WORKFLOW_ID=4

# Update dataset info
poetry run client dataset modify-dataset $PROJECT_ID $DATASET_ID --name_dataset test-ds --type image --read_only true

# Add resource to dataset
poetry run client dataset add-resource $PROJECT_ID $DATASET_ID ${TESTDATA}/png/ --glob_pattern *.png

# Add output dataset
poetry run client project add-dataset $PROJECT_ID out-ds --type zarr
poetry run client dataset add-resource $PROJECT_ID $DATASET_ID_OUT tmp-proj/out-ds --glob_pattern *.zarr


# These two lines are taken care of by alembic upgrade
#poetry run client task new "Create OME-ZARR structure" task img zarr
#poetry run client task new "Yokogawa to Zarr" zarr zarr


# Create workflow
poetry run client task new "My WF" workflow image zarr

# Add subtasks (with args, if needed)
poetry run client task add-subtask "My WF" "Create OME-ZARR structure"
echo "{\"parallelization_level\" : \"well\", \"rows\":2, \"cols\": 1}" > tmp-proj/args_yoko.json
poetry run client task add-subtask "My WF" "Yokogawa to Zarr" --args_json tmp-proj/args_yoko.json

# Apply workflow
poetry run client workflow apply $PROJECT_ID $DATASET_ID $WORKFLOW_ID --output_dataset_id $DATASET_ID_OUT

# Show info
poetry run client project list
poetry run client dataset show $PROJECT_ID test-ds
poetry run client dataset show-resources $PROJECT_ID $DATASET_ID
poetry run client dataset show $PROJECT_ID out-ds
poetry run client dataset show-resources $PROJECT_ID $DATASET_ID_OUT
poetry run client task list
