LABEL=2
TMPDIR=tmp-proj-$LABEL

rm -r $TMPDIR
mkdir $TMPDIR

# Create project
poetry run client project new test-prj $TMPDIR

TESTDATA=../../tests/data
WFNAME="My WF $LABEL"


PROJECT_ID=2
DATASET_ID=2
DATASET_ID_OUT=4
WORKFLOW_ID=4

# Update dataset info
poetry run client dataset modify-dataset $PROJECT_ID $DATASET_ID --name_dataset test-ds-$LABEL --type image --read_only true

# Add resource to dataset
poetry run client dataset add-resource $PROJECT_ID $DATASET_ID ${TESTDATA}/png/ --glob_pattern *.png

# Add output dataset
poetry run client project add-dataset $PROJECT_ID out-ds-$LABEL --type zarr
poetry run client dataset add-resource $PROJECT_ID $DATASET_ID_OUT ${TMPDIR}/out-ds-$LABEL --glob_pattern *.zarr


# Create workflow
poetry run client task new "$WFNAME" workflow image zarr

# Add subtasks (with args, if needed)
poetry run client task add-subtask "$WFNAME" "Create OME-ZARR structure"
echo "{\"parallelization_level\" : \"well\", \"rows\":2, \"cols\": 1}" > ${TMPDIR}/args_yoko.json
poetry run client task add-subtask "$WFNAME" "Yokogawa to Zarr" --args_json ${TMPDIR}/args_yoko.json

# Apply workflow
poetry run client workflow apply $PROJECT_ID $DATASET_ID $WORKFLOW_ID --output_dataset_id $DATASET_ID_OUT

# Show info
poetry run client project list
poetry run client dataset show $PROJECT_ID test-ds-$LABEL
poetry run client dataset show-resources $PROJECT_ID $DATASET_ID
poetry run client dataset show $PROJECT_ID out-ds-$LABEL
poetry run client dataset show-resources $PROJECT_ID $DATASET_ID_OUT
poetry run client task list
