# Register user (this step will change in the future)
http POST localhost:8000/auth/register email=test@me.com password=test

# Define/initialize empty folder for project-related info
# (and also for the output dataset -- see below)
TMPDIR=`pwd`/tmp-proj-1
rm -r $TMPDIR
mkdir $TMPDIR

# Set useful variables
PROJECT_NAME="project_10x10"
DATASET_IN_NAME="input-ds-1"
DATASET_OUT_NAME="output-ds-1"
WORKFLOW_NAME="My worfklow 1"

# Create project
fractal project new $PROJECT_NAME $TMPDIR

INPUT_PATH=/data/active/fractal/3D/PelkmansLab/CardiacMultiplexing/Cycle1_5x5_10wells_constantZ
OUTPUT_PATH=/data/active/fractal/tests/20220919_10well_5x5

# Update dataset info
fractal dataset modify-dataset $PROJECT_NAME "default" --new_dataset_name $DATASET_IN_NAME --type image --read_only true

# Add resource to dataset
fractal dataset add-resource $PROJECT_NAME $DATASET_IN_NAME ${INPUT_PATH} --glob_pattern *.png

# Add output dataset
fractal project add-dataset $PROJECT_NAME $DATASET_OUT_NAME --type zarr
fractal dataset add-resource $PROJECT_NAME $DATASET_OUT_NAME ${OUTPUT_PATH} --glob_pattern *.zarr


# Create workflow
fractal task new "$WORKFLOW_NAME" workflow image zarr

# Add subtasks (with args, if needed)
echo "{\"num_levels\": 5, \"coarsening_xy\": 2, \"channel_parameters\": {\"A01_C01\": {\"label\": \"DAPI\",\"colormap\": \"00FFFF\",\"start\": 110,\"end\": 800 }, \"A01_C02\": {\"label\": \"nanog\",\"colormap\": \"FF00FF\",\"start\": 110,\"end\": 290 }, \"A02_C03\": {\"label\": \"Lamin B1\",\"colormap\": \"FFFF00\",\"start\": 110,\"end\": 1600 }}}" > ${TMPDIR}/args_create.json
fractal task add-subtask "$WORKFLOW_NAME" "Create OME-ZARR structure" --args_json ${TMPDIR}/args_create.json

echo "{\"executor\": \"cpu\"}" > ${TMPDIR}/args_yoko.json
fractal task add-subtask "$WORKFLOW_NAME" "Yokogawa to Zarr" --args_json ${TMPDIR}/args_yoko.json

cp params_uzh_cardiac_illumination.json ${TMPDIR}/args_illum.json
fractal task add-subtask "$WORKFLOW_NAME" "Illumination correction" --args_json ${TMPDIR}/args_illum.json

echo "{\"labeling_level\": 1, \"labeling_channel\": \"A01_C01\", \"executor\": \"gpu\"}" > ${TMPDIR}/args_labeling.json
fractal task add-subtask "$WORKFLOW_NAME" "Per-FOV image labeling" --args_json ${TMPDIR}/args_labeling.json

fractal task add-subtask "$WORKFLOW_NAME" "Replicate Zarr structure"
echo "{\"executor\": \"cpu\"}" > ${TMPDIR}/args_mip.json
fractal task add-subtask "$WORKFLOW_NAME" "Maximum Intensity Projection" --args_json ${TMPDIR}/args_mip.json

echo "{\"labeling_level\": 2, \"labeling_channel\": \"A01_C01\", \"executor\": \"gpu\"}" > ${TMPDIR}/args_whole_well_labeling.json
fractal task add-subtask "$WORKFLOW_NAME" "Whole-well image labeling" --args_json ${TMPDIR}/args_whole_well_labeling.json

echo "{\"level\": 0, \"table_name\": \"nuclei\", \"executor\": \"cpu\", \"workflow_file\": \"/data/homes/jluethi/fractal_3repo/fractal/examples/05_10x10_test_constant_z/regionprops_from_existing_labels_feature.yaml\"}" > ${TMPDIR}/args_measurement.json
fractal task add-subtask "$WORKFLOW_NAME" "Measurement" --args_json ${TMPDIR}/args_measurement.json

# Apply workflow
fractal workflow apply $PROJECT_NAME $DATASET_IN_NAME "$WORKFLOW_NAME" --output_dataset_name $DATASET_OUT_NAME
