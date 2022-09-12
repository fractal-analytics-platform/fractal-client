# Register user (this step will change in the future)
http POST localhost:8000/auth/register email=test@me.com password=test

# Define/initialize empty folder for project-related info
# (and also for the output dataset -- see below)
TMPDIR=`pwd`/tmp-proj-1
rm -r $TMPDIR
mkdir $TMPDIR

# Set useful variables
PROJECT_NAME="myproj-1"
DATASET_IN_NAME="input-ds-1"
DATASET_OUT_NAME="output-ds-1"
WORKFLOW_NAME="My worfklow 1"

# Create project
fractal project new $PROJECT_NAME $TMPDIR

INPUT_PATH=../images/10.5281_zenodo.7059515

# Update dataset info
fractal dataset modify-dataset $PROJECT_NAME "default" --new_dataset_name $DATASET_IN_NAME --type image --read_only true

# Add resource to dataset
fractal dataset add-resource $PROJECT_NAME $DATASET_IN_NAME ${INPUT_PATH} --glob_pattern *.png

# Add output dataset
fractal project add-dataset $PROJECT_NAME $DATASET_OUT_NAME --type zarr
fractal dataset add-resource $PROJECT_NAME $DATASET_OUT_NAME ${TMPDIR}/$DATASET_OUT_NAME --glob_pattern *.zarr


# Create workflow
fractal task new "$WORKFLOW_NAME" workflow image zarr

# Add subtasks (with args, if needed)
echo "{\"num_levels\": 5, \"coarsening_xy\": 2}" > ${TMPDIR}/args_create.json
fractal task add-subtask "$WORKFLOW_NAME" "Create OME-ZARR structure" --args_json ${TMPDIR}/args_create.json

echo "{\"parallelization_level\" : \"well\", \"rows\":1, \"cols\": 2, \"executor\": \"cpu\"}" > ${TMPDIR}/args_yoko.json
fractal task add-subtask "$WORKFLOW_NAME" "Yokogawa to Zarr" --args_json ${TMPDIR}/args_yoko.json

echo "{\"parallelization_level\" : \"well\", \"labeling_level\": 1, \"executor\": \"gpu\"}" > ${TMPDIR}/args_labeling.json
fractal task add-subtask "$WORKFLOW_NAME" "Per-FOV image labeling" --args_json ${TMPDIR}/args_labeling.json

fractal task add-subtask "$WORKFLOW_NAME" "Replicate Zarr structure"
echo "{\"parallelization_level\" : \"well\", \"executor\": \"cpu\"}" > ${TMPDIR}/args_mip.json
fractal task add-subtask "$WORKFLOW_NAME" "Maximum Intensity Projection" --args_json ${TMPDIR}/args_mip.json

# Apply workflow
fractal workflow apply $PROJECT_NAME $DATASET_IN_NAME "$WORKFLOW_NAME" --output_dataset_name $DATASET_OUT_NAME
