# Register user (this step will change in the future)
curl -d '{"email":"test@me.com", "password":"test"}' -H "Content-Type: application/json" -X POST localhost:8000/auth/register

# Set useful variables
PRJ_NAME="myproj-1w-2x2-3c-10z"
DS_IN_NAME="input-ds"
DS_OUT_NAME="output-ds"
WF_NAME="WF 1w-2x2-3c-10z"
export FRACTAL_CACHE_PATH=`pwd`/".cache"
rm -v ${FRACTAL_CACHE_PATH}/session
rm -v ${FRACTAL_CACHE_PATH}/tasks

# Define/initialize empty folder for temporary files
TMPDIR=`pwd`/$PRJ_NAME
rm -r $TMPDIR
mkdir $TMPDIR
TMPJSON=${TMPDIR}/tmp.json

INPUT_PATH=../images/10.5281_zenodo.7057076
OUTPUT_PATH=${TMPDIR}/output

CMD="fractal"
CMD_JSON="python aux_extract_id_from_project_json.py $TMPJSON"

# Create project
$CMD -j project new $PRJ_NAME $TMPDIR > $TMPJSON
PRJ_ID=`$CMD_JSON project_id`
DS_IN_ID=`$CMD_JSON dataset_id "default"`
echo "PRJ_ID: $PRJ_ID"
echo "DS_IN_ID: $DS_IN_ID"

# Update dataset name/type, and add a resource
$CMD dataset edit --name "$DS_IN_NAME" -t image --read-only $PRJ_ID $DS_IN_ID
$CMD dataset add-resource -g "*.png" $PRJ_ID $DS_IN_ID $INPUT_PATH

# Add output dataset, and add a resource to it
DS_OUT_ID=`$CMD --batch project add-dataset $PRJ_ID "$DS_OUT_NAME"`
$CMD dataset edit -t zarr --read-write $PRJ_ID $DS_OUT_ID
$CMD dataset add-resource -g "*.zarr" $PRJ_ID $DS_OUT_ID $OUTPUT_PATH

# Create workflow
WF_ID=`$CMD --batch task new "$WF_NAME" workflow image zarr`
echo "WF_ID: $WF_ID"

# Add subtasks

echo "{\"num_levels\": 5, \"coarsening_xy\": 2, \"channel_parameters\": {\"A01_C01\": {\"label\": \"DAPI\",\"colormap\": \"00FFFF\",\"start\": 110,\"end\": 800 }, \"A01_C02\": {\"label\": \"nanog\",\"colormap\": \"FF00FF\",\"start\": 110,\"end\": 290 }, \"A02_C03\": {\"label\": \"Lamin B1\",\"colormap\": \"FFFF00\",\"start\": 110,\"end\": 1600 }}}" > ${TMPDIR}/args_create.json
$CMD task add-subtask $WF_ID "Create OME-ZARR structure" --args-file ${TMPDIR}/args_create.json

$CMD task add-subtask $WF_ID "Yokogawa to Zarr"

echo "{\"labeling_level\": 3, \"executor\": \"gpu\"}" > ${TMPDIR}/args_labeling.json
$CMD task add-subtask $WF_ID "Cellpose Segmentation" --args-file ${TMPDIR}/args_labeling.json

$CMD task add-subtask $WF_ID "Replicate Zarr structure"

$CMD task add-subtask $WF_ID "Maximum Intensity Projection"

# Apply workflow
$CMD task apply $PRJ_ID $DS_IN_ID $DS_OUT_ID $WF_ID
