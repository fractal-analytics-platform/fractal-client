# Register user (this step will change in the future)
http POST localhost:8000/auth/register email=test@me.com password=test

# Define/initialize empty folder for project-related info
# (and also for the output dataset -- see below)
CURRENT_DIR=`pwd`
TMPDIR=${CURRENT_DIR}/tmp-proj-ww
rm -r $TMPDIR
mkdir $TMPDIR

# Set useful variables
PROJECT_NAME="myproj-ww"
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

echo "{\"__PROVIDER_ARGS__\" : {\"max_blocks\": 10}}" > /tmp/args_wf.json
poetry run client task modify-task "$WORKFLOW_NAME" --default_args /tmp/args_wf.json

# Add subtasks (with args, if needed)
echo "{\"num_levels\": 4}" > /tmp/args_create.json
poetry run client task add-subtask "$WORKFLOW_NAME" "Create OME-ZARR structure" --args_json /tmp/args_create.json

echo "{\"parallelization_level\" : \"well\", \"rows\":1, \"cols\": 2}" > /tmp/args_yoko.json
poetry run client task add-subtask "$WORKFLOW_NAME" "Yokogawa to Zarr" --args_json /tmp/args_yoko.json

poetry run client task add-subtask "$WORKFLOW_NAME" "Replicate Zarr structure"

echo "{\"parallelization_level\" : \"well\"}" > /tmp/args_mip.json
poetry run client task add-subtask "$WORKFLOW_NAME" "Maximum Intensity Projection" --args_json /tmp/args_mip.json

echo "{\"parallelization_level\" : \"well\", \"labeling_level\": 3}" > /tmp/args_labeling.json
poetry run client task add-subtask "$WORKFLOW_NAME" "Whole-well image labeling" --args_json /tmp/args_labeling.json

echo "{\"parallelization_level\" : \"well\",\
       \"workflow_file\": \"${CURRENT_DIR}/regionprops_from_existing_labels.yaml\",\
       \"table_name\": \"nuclei\",\
       \"whole_well\": true}" > /tmp/args_measurement.json
poetry run client task add-subtask "$WORKFLOW_NAME" "Measurement" --args_json /tmp/args_measurement.json

# Apply workflow
poetry run client workflow apply $PROJECT_NAME $DATASET_IN_NAME "$WORKFLOW_NAME" --output_dataset_name $DATASET_OUT_NAME
