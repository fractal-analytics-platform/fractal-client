# Register user (this step will change in the future)
curl -d '{"email":"test@me.com", "password":"test"}' -H "Content-Type: application/json" -X POST localhost:8000/auth/register
echo

# Set useful variables
PRJ_NAME="myproj-zenodo-2x2_mip"
DS_IN_NAME="input-ds2"
DS_OUT_NAME="output-ds2"
WF_NAME="WF zenodo 2x2_mip"

# Define/initialize empty folder for temporary files
TMPDIR=`pwd`/$PRJ_NAME
rm -r $TMPDIR
mkdir $TMPDIR

INPUT_PATH=../images/10.5281_zenodo.7057076
OUTPUT_PATH=/data/active/jluethi/Fractal/20220928_2x2_Zenodo
rm -rv $OUTPUT_PATH

TMPJSON=${TMPDIR}/tmp.json
TMPTASKS=${TMPDIR}/core_tasks.json

CMD="fractal"
CMD_JSON="python aux_extract_from_simple_json.py $TMPJSON"
CMD_CORE_TASKS="python aux_extract_id_for_core_task.py $TMPTASKS"
fractal task list > $TMPTASKS

# Create project
fractal -j project new $PRJ_NAME $TMPDIR > $TMPJSON
PRJ_ID=`$CMD_JSON id`
DS_IN_ID=`$CMD_JSON id`
echo "PRJ_ID: $PRJ_ID"
echo "DS_IN_ID: $DS_IN_ID"

# Update dataset name/type, and add a resource
fractal dataset edit --name "$DS_IN_NAME" -t image --read-only $PRJ_ID $DS_IN_ID
fractal dataset add-resource -g "*.png" $PRJ_ID $DS_IN_ID $INPUT_PATH

# Add output dataset, and add a resource to it
DS_OUT_ID=`fractal --batch project add-dataset $PRJ_ID "$DS_OUT_NAME"`
fractal dataset edit -t zarr --read-write $PRJ_ID $DS_OUT_ID
fractal dataset add-resource -g "*.zarr" $PRJ_ID $DS_OUT_ID $OUTPUT_PATH

# Create workflow
WF_ID=`fractal --batch task new "$WF_NAME" workflow image zarr`
echo "WF_ID: $WF_ID"

# Add subtasks

SUBTASK_ID=`$CMD_CORE_TASKS "Create OME-ZARR structure"`
echo "{\"num_levels\": 5, \"coarsening_xy\": 2, \"channel_parameters\": {\"A01_C01\": {\"label\": \"DAPI\",\"colormap\": \"00FFFF\",\"start\": 50,\"end\": 700 }, \"A01_C02\": {\"label\": \"nanog\",\"colormap\": \"FF00FF\",\"start\": 20,\"end\": 200 }, \"A02_C03\": {\"label\": \"Lamin B1\",\"colormap\": \"FFFF00\",\"start\": 50,\"end\": 1500 }}}" > ${TMPDIR}/args_create.json
fractal task add-subtask $WF_ID $SUBTASK_ID --args-file ${TMPDIR}/args_create.json

SUBTASK_ID=`$CMD_CORE_TASKS "Yokogawa to Zarr"`
fractal task add-subtask $WF_ID $SUBTASK_ID

SUBTASK_ID=`$CMD_CORE_TASKS "Illumination correction"`
# Paths of illumination correction images need to be accessible on the server.
# This works if one runs the client from the same machine as the server. Otherwise, change `root_path_corr`
echo "{\"overwrite\": true, \"executor\": \"cpu-mid\", \"dict_corr\": {\"root_path_corr\": \"$TMPDIR/../illum_corr_images/\", \"A01_C01\": \"20220621_UZH_manual_illumcorr_40x_A01_C01.png\", \"A01_C02\": \"20220621_UZH_manual_illumcorr_40x_A01_C02.png\", \"A02_C03\": \"20220621_UZH_manual_illumcorr_40x_A02_C03.png\"}}" > ${TMPDIR}/args_illum.json
fractal task add-subtask $WF_ID $SUBTASK_ID --args-file ${TMPDIR}/args_illum.json

SUBTASK_ID=`$CMD_CORE_TASKS "Replicate Zarr structure"`
echo "{\"executor\": \"cpu-mid\"}" > ${TMPDIR}/args_replicate.json
fractal task add-subtask $WF_ID $SUBTASK_ID --args-file ${TMPDIR}/args_replicate.json

SUBTASK_ID=`$CMD_CORE_TASKS "Maximum Intensity Projection"`
echo "{\"executor\": \"cpu-mid\"}" > ${TMPDIR}/args_mip.json
fractal task add-subtask $WF_ID $SUBTASK_ID --args-file ${TMPDIR}/args_mip.json

SUBTASK_ID=`$CMD_CORE_TASKS "Cellpose Segmentation"`
echo "{\"labeling_level\": 2, \"executor\": \"cpu-mid\", \"ROI_table_name\": \"well_ROI_table\"}" > ${TMPDIR}/args_labeling.json
fractal task add-subtask $WF_ID $SUBTASK_ID --args-file ${TMPDIR}/args_labeling.json

SUBTASK_ID=`$CMD_CORE_TASKS "Measurement"`
echo "{\"level\": 0, \"measurement_table_name\": \"nuclei\", \"executor\": \"cpu-mid\", \"ROI_table_name\": \"well_ROI_table\",\"workflow_file\": \"$TMPDIR/../regionprops_from_existing_labels_feature.yaml\"}" > ${TMPDIR}/args_measurement.json
fractal task add-subtask $WF_ID $SUBTASK_ID --args-file ${TMPDIR}/args_measurement.json

# Apply workflow
fractal task apply $PRJ_ID $DS_IN_ID $DS_OUT_ID $WF_ID
