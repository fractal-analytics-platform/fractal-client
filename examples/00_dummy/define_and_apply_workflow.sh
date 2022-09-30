# Register user (this step will change in the future)
curl -d '{"email":"test@me.com", "password":"test"}' -H "Content-Type: application/json" -X POST localhost:8000/auth/register

# Set useful variables
PRJ_NAME="myproj-dummy"
DS_IN_NAME="input-ds-dummy"
DS_OUT_NAME="output-ds-dummy"
WF_NAME="My workflow dummy"
export FRACTAL_CACHE_PATH=`pwd`/".cache"

# Define/initialize empty folder for temporary files
TMPDIR=`pwd`/$PRJ_NAME
rm -r $TMPDIR
mkdir $TMPDIR
TMPJSON=${TMPDIR}/tmp.json

INPUT_PATH=/tmp
OUTPUT_PATH=${TMPDIR}/output

CMD="fractal"
CMD_JSON="python aux_extract_from_simple_json.py $TMPJSON"

# Create project
$CMD -j project new $PRJ_NAME $TMPDIR > $TMPJSON
PRJ_ID=`$CMD_JSON id`
DS_IN_ID=`$CMD_JSON id`
echo "PRJ_ID: $PRJ_ID"
echo "DS_IN_ID: $DS_IN_ID"

# Update dataset name/type, and add a resource
$CMD dataset edit --name "$DS_IN_NAME" -t image --read-only $PRJ_ID $DS_IN_ID
$CMD dataset add-resource -g "*.png" $PRJ_ID $DS_IN_ID $INPUT_PATH

# Add output dataset, and add a resource to it
DS_OUT_ID=`$CMD --batch project add-dataset $PRJ_ID "$DS_OUT_NAME"`
$CMD dataset edit -t none --read-write $PRJ_ID $DS_OUT_ID
$CMD dataset add-resource -g "*.json" $PRJ_ID $DS_OUT_ID $OUTPUT_PATH

# Create workflow
WF_ID=`$CMD --batch task new "$WF_NAME" workflow image zarr`
echo "WF_ID: $WF_ID"

# Add subtasks
echo "{\"message\": \"bottle\", \"index\": 1}" > ${TMPDIR}/args_tmp.json
$CMD task add-subtask $WF_ID "dummy" --args-file ${TMPDIR}/args_tmp.json

# Apply workflow
$CMD task apply $PRJ_ID $DS_IN_ID $DS_OUT_ID $WF_ID
