# Register user (this step will change in the future)
http POST localhost:8000/auth/register email=test@me.com password=test

CMD="poetry run fractal"

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
$CMD project new $PROJECT_NAME $TMPDIR

INPUT_PATH=../images/10.5281_zenodo.7059515

# Update dataset info
$CMD dataset modify-dataset $PROJECT_NAME "default" --new_dataset_name $DATASET_IN_NAME --type image --read_only true

# Add resource to dataset
$CMD dataset add-resource $PROJECT_NAME $DATASET_IN_NAME ${INPUT_PATH} --glob_pattern *.png

# Add output dataset
$CMD project add-dataset $PROJECT_NAME $DATASET_OUT_NAME --type zarr
$CMD dataset add-resource $PROJECT_NAME $DATASET_OUT_NAME ${TMPDIR}/$DATASET_OUT_NAME --glob_pattern *.zarr


# Create workflow
$CMD task new "$WORKFLOW_NAME" workflow image zarr

echo "{}" > /tmp/args_wf.json
$CMD task modify-task "$WORKFLOW_NAME" --default_args /tmp/args_wf.json

# Add subtasks (with args, if needed)
echo "{\"num_levels\": 4, \"coarsening_xy\": 2}" > /tmp/args_create.json
$CMD task add-subtask "$WORKFLOW_NAME" "Create OME-ZARR structure" --args_json /tmp/args_create.json

echo "{\"parallelization_level\" : \"well\", \"rows\":1, \"cols\": 2}" > /tmp/args_yoko.json
$CMD task add-subtask "$WORKFLOW_NAME" "Yokogawa to Zarr" --args_json /tmp/args_yoko.json

echo "{\"parallelization_level\" : \"well\", \"labeling_level\": 3}" > /tmp/args_labeling.json
$CMD task add-subtask "$WORKFLOW_NAME" "Per-FOV image labeling" --args_json /tmp/args_labeling.json

$CMD task add-subtask "$WORKFLOW_NAME" "Replicate Zarr structure"
echo "{\"parallelization_level\" : \"well\"}" > /tmp/args_mip.json
$CMD task add-subtask "$WORKFLOW_NAME" "Maximum Intensity Projection" --args_json /tmp/args_mip.json

# Apply workflow
$CMD workflow apply $PROJECT_NAME $DATASET_IN_NAME "$WORKFLOW_NAME" --output_dataset_name $DATASET_OUT_NAME
