http POST localhost:8000/auth/register email=test2@me.com password=test2


LABEL=2
TMPDIR=user$LABEL/tmp-proj-$LABEL

rm -r ../$TMPDIR
mkdir ../$TMPDIR

PROJECT_NAME="myproj"
DATASET_IN_NAME="in-ds-$LABEL"
DATASET_OUT_NAME="out-ds-$LABEL"
WORKFLOW_NAME="My WF $LABEL"

# Create project
poetry run client project new $PROJECT_NAME $TMPDIR

TESTDATA=../../tests/data

# Update dataset info
poetry run client dataset modify-dataset $PROJECT_NAME "default" --new_dataset_name $DATASET_IN_NAME --type image --read_only true

# Add resource to dataset
poetry run client dataset add-resource $PROJECT_NAME $DATASET_IN_NAME ${TESTDATA}/png/ --glob_pattern *.png

# Add output dataset
poetry run client project add-dataset $PROJECT_NAME $DATASET_OUT_NAME --type zarr
poetry run client dataset add-resource $PROJECT_NAME $DATASET_OUT_NAME ${TMPDIR}/$DATASET_OUT_NAME --glob_pattern *.zarr


# Create workflow
poetry run client task new "$WORKFLOW_NAME" workflow image zarr

echo "{\"__PROVIDER_ARGS__\" : {\"max_blocks\": 10}}" > /tmp/args_wf_${LABEL}.json
poetry run client task modify-task "$WORKFLOW_NAME" --default_args /tmp/args_wf_${LABEL}.json

# Add subtasks (with args, if needed)
poetry run client task add-subtask "$WORKFLOW_NAME" "Create OME-ZARR structure"
echo "{\"parallelization_level\" : \"well\", \"rows\":2, \"cols\": 1}" > /tmp/args_yoko_${LABEL}.json
poetry run client task add-subtask "$WORKFLOW_NAME" "Yokogawa to Zarr" --args_json /tmp/args_yoko_${LABEL}.json

# Apply workflow
poetry run client workflow apply $PROJECT_NAME $DATASET_IN_NAME "$WORKFLOW_NAME" --output_dataset_name $DATASET_OUT_NAME

# Show info
poetry run client project list
poetry run client dataset show $PROJECT_NAME $DATASET_IN_NAME
poetry run client dataset show-resources $PROJECT_NAME $DATASET_IN_NAME
poetry run client dataset show $PROJECT_NAME $DATASET_OUT_NAME
poetry run client dataset show-resources $PROJECT_NAME $DATASET_OUT_NAME
poetry run client task list
