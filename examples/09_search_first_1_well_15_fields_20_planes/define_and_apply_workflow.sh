# Register user (this step will change in the future)
curl -d '{"email":"test@me.com", "password":"test"}' -H "Content-Type: application/json" -X POST localhost:8000/auth/register
echo

# Set useful variables
PRJ_NAME="myproj-10w-5x5"
DS_IN_NAME="input-ds"
DS_OUT_NAME="output-ds"
WF_NAME="WF 10w-5x5"

# Define/initialize empty folder for temporary files
TMPDIR=`pwd`/$PRJ_NAME
rm -r $TMPDIR
mkdir $TMPDIR

INPUT_PATH=/data/active/fractal/Liberali/1_well_15_fields_20_planes_SF_w_errors/D10_R1/220304_172545_220304_175557
OUTPUT_PATH=/data/active/jluethi/Fractal/20220926_SearchFirst_15
rm -rv $OUTPUT_PATH

TMPJSON=${TMPDIR}/tmp.json
TMPTASKS=${TMPDIR}/core_tasks.json

CMD="fractal"
CMD_JSON="python aux_extract_from_simple_json.py $TMPJSON"
CMD_CORE_TASKS="python aux_extract_id_for_core_task.py $TMPTASKS"
$CMD task list > $TMPTASKS

# Create project
$CMD -j project new $PRJ_NAME $TMPDIR > $TMPJSON
PRJ_ID=`$CMD_JSON id`
DS_IN_ID=`$CMD_JSON id`
echo "PRJ_ID: $PRJ_ID"
echo "DS_IN_ID: $DS_IN_ID"

# Update dataset name/type, and add a resource
$CMD dataset edit --name "$DS_IN_NAME" -t image --read-only $PRJ_ID $DS_IN_ID
$CMD dataset add-resource -g "*.tif" $PRJ_ID $DS_IN_ID $INPUT_PATH

# Add output dataset, and add a resource to it
DS_OUT_ID=`$CMD --batch project add-dataset $PRJ_ID "$DS_OUT_NAME"`
$CMD dataset edit -t zarr --read-write $PRJ_ID $DS_OUT_ID
$CMD dataset add-resource -g "*.zarr" $PRJ_ID $DS_OUT_ID $OUTPUT_PATH

# Create workflow
WF_ID=`$CMD --batch task new "$WF_NAME" workflow image zarr`
echo "WF_ID: $WF_ID"

# Add subtasks

SUBTASK_ID=`$CMD_CORE_TASKS "Create OME-ZARR structure"`
echo "{\"num_levels\": 5, \"coarsening_xy\": 2, \"channel_parameters\": {\"A01_C01\": {\"label\": \"Channel 1\",\"colormap\": \"00FFFF\",\"start\": 110,\"end\": 2000 }, \"A02_C02\": {\"label\": \"Channel 2\",\"colormap\": \"FF00FF\",\"start\": 110,\"end\": 500}, \"A03_C03\": {\"label\": \"Channel 3\",\"colormap\": \"00FF00\",\"start\": 110,\"end\": 1600 }, \"A04_C04\": {\"label\": \"Channel 4\",\"colormap\": \"FFFF00\",\"start\": 110,\"end\": 1600 }}}" > ${TMPDIR}/args_create.json
$CMD task add-subtask $WF_ID $SUBTASK_ID --args-file ${TMPDIR}/args_create.json

SUBTASK_ID=`$CMD_CORE_TASKS "Yokogawa to Zarr"`
$CMD task add-subtask $WF_ID $SUBTASK_ID

# SUBTASK_ID=`$CMD_CORE_TASKS "Illumination correction"`
# echo "{\"overwrite\": true, \"executor\": \"cpu-mid\", \"dict_corr\": {\"root_path_corr\": \"/data/active/fractal/3D/PelkmansLab/IlluminationCorrection_Matrices_UZH/\", \"A01_C01\": \"20220621_UZH_manual_illumcorr_40x_A01_C01.tif\", \"A01_C02\": \"20220621_UZH_manual_illumcorr_40x_A01_C02.tif\", \"A02_C03\": \"20220621_UZH_manual_illumcorr_40x_A02_C03.tif\"}}" > ${TMPDIR}/args_illum.json
# $CMD task add-subtask $WF_ID $SUBTASK_ID --args-file ${TMPDIR}/args_illum.json

SUBTASK_ID=`$CMD_CORE_TASKS "Replicate Zarr structure"`
echo "{\"executor\": \"cpu-low\"}" > ${TMPDIR}/args_replicate.json
$CMD task add-subtask $WF_ID $SUBTASK_ID --args-file ${TMPDIR}/args_replicate.json

SUBTASK_ID=`$CMD_CORE_TASKS "Maximum Intensity Projection"`
echo "{\"executor\": \"cpu-mid\"}" > ${TMPDIR}/args_mip.json
$CMD task add-subtask $WF_ID $SUBTASK_ID --args-file ${TMPDIR}/args_mip.json

SUBTASK_ID=`$CMD_CORE_TASKS "Cellpose Segmentation"`
echo "{\"labeling_level\": 3, \"executor\": \"cpu-mid\", \"ROI_table_name\": \"well_ROI_table\", \"diameter_level0\": 1200.0, \"cellprob_threshold\": -1.0, \"flow_threshold\":  0.6}" > ${TMPDIR}/args_labeling.json
$CMD task add-subtask $WF_ID $SUBTASK_ID --args-file ${TMPDIR}/args_labeling.json

SUBTASK_ID=`$CMD_CORE_TASKS "Measurement"`
echo "{\"level\": 0, \"measurement_table_name\": \"organoid_measurements\", \"ROI_table_name\": \"well_ROI_table\", \"executor\": \"cpu-mid\", \"workflow_file\": \"${TMPDIR}/../regionprops_from_existing_labels_feature.yaml\"}" > ${TMPDIR}/args_measurement.json
$CMD task add-subtask $WF_ID $SUBTASK_ID --args-file ${TMPDIR}/args_measurement.json

# Apply workflow
$CMD task apply $PRJ_ID $DS_IN_ID $DS_OUT_ID $WF_ID
