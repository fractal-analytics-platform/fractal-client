# Register user (this step will change in the future)
http POST localhost:8000/auth/register email=test@me.com password=test

# Define/initialize empty folder for project-related info
# (and also for the output dataset -- see below)
TMPDIR=`pwd`/tmp-proj
rm -r $TMPDIR
mkdir $TMPDIR

# Set useful variables
PROJECT_NAME="myproj"
DATASET_IN_NAME="input-ds"
DATASET_OUT_NAME="output-ds"
WORKFLOW_NAME="My worfklow"

# Create project
poetry run client project new $PROJECT_NAME $TMPDIR

TESTDATA=../../../tests/data

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

echo "{\"parallelization_level\" : \"well\", \"rows\":1, \"cols\": 2}" > /tmp/args_yoko_${LABEL}.json
poetry run client task add-subtask "$WORKFLOW_NAME" "Yokogawa to Zarr" --args_json /tmp/args_yoko_${LABEL}.json

poetry run client task add-subtask "$WORKFLOW_NAME" "Replicate Zarr structure"

echo "{\"parallelization_level\" : \"well\"}" > /tmp/args_mip_${LABEL}.json
poetry run client task add-subtask "$WORKFLOW_NAME" "Maximum Intensity Projection" --args_json /tmp/args_mip_${LABEL}.json

echo "{\"parallelization_level\" : \"well\"}" > /tmp/args_labeling_${LABEL}.json
poetry run client task add-subtask "$WORKFLOW_NAME" "Per-FOV image labeling" --args_json /tmp/args_labeling_${LABEL}.json

# Apply workflow
poetry run client workflow apply $PROJECT_NAME $DATASET_IN_NAME "$WORKFLOW_NAME" --output_dataset_name $DATASET_OUT_NAME
