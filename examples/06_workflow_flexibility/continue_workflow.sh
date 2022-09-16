# Set useful variables
TMPDIR=`pwd`/tmp-proj-1

PROJECT_NAME="myproj-1"
DATASET_OUT_NAME="output-ds-1"
WORKFLOW_NAME="Continuation Workflow"

# Create workflow
fractal task new "$WORKFLOW_NAME" workflow zarr zarr

fractal task add-subtask "$WORKFLOW_NAME" "Replicate Zarr structure"
echo "{\"parallelization_level\" : \"well\", \"executor\": \"cpu\"}" > ${TMPDIR}/args_mip.json
fractal task add-subtask "$WORKFLOW_NAME" "Maximum Intensity Projection" --args_json ${TMPDIR}/args_mip.json

echo "{\"parallelization_level\" : \"well\", \"labeling_level\": 3, \"labeling_channel\": \"A01_C01\", \"executor\": \"gpu\"}" > ${TMPDIR}/args_whole_well_labeling.json
fractal task add-subtask "$WORKFLOW_NAME" "Whole-well image labeling" --args_json ${TMPDIR}/args_whole_well_labeling.json

# TODO: Couldn't make relative path to regionprops_from_existing_labels_feature.yaml work, currently it's the absolute path
#echo "{\"parallelization_level\" : \"well\", \"level\": 0, \"table_name\": \"nuclei\", \"executor\": \"cpu\", \"workflow_file\": \"/data/homes/jluethi/fractal_3repo/fractal/examples/06_workflow_flexibility/regionprops_from_existing_labels_feature.yaml\"}" > ${TMPDIR}/args_measurement.json
echo "{\"parallelization_level\" : \"well\", \"level\": 0, \"table_name\": \"nuclei\", \"executor\": \"cpu\", \"workflow_file\": \"$TMPDIR/../regionprops_from_existing_labels_feature.yaml\"}" > ${TMPDIR}/args_measurement.json
fractal task add-subtask "$WORKFLOW_NAME" "Measurement" --args_json ${TMPDIR}/args_measurement.json

# Apply workflow
fractal workflow apply $PROJECT_NAME $DATASET_OUT_NAME "$WORKFLOW_NAME" --output_dataset_name $DATASET_OUT_NAME
