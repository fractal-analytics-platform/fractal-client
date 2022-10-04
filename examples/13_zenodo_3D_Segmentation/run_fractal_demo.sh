# Register user (this step will change in the future)
curl -d '{"email":"test@me.com", "password":"test"}' -H "Content-Type: application/json" -X POST localhost:8000/auth/register
echo -e "FRACTAL_USER=test@me.com\nFRACTAL_PASSWORD=test" > .fractal.env
echo

# Temporary: Set a unique name prefix to ensure all names are unique
# This will not be necessary long-term, but is a workaround for some current constraints
# Plus we use it in the demo to define unique output folders
USERNAME='joel'
if [ "$USERNAME" = "" ];
then
    echo "Please define a username on line 8"
    exit 1
fi

# Set useful variables
PRJ_NAME=$USERNAME"_myproj"
DS_IN_NAME=$USERNAME"_input-ds"
DS_OUT_NAME=$USERNAME"_output-ds"
WF_NAME=$USERNAME"_WF-2x2"
export FRACTAL_CACHE_PATH=`pwd`/".cache"

# Define/initialize empty folder for temporary files
TMPDIR=`pwd`/$PRJ_NAME
rm -r $TMPDIR
mkdir $TMPDIR

# If the images have not been downloaded yet, use the `fetch_test_data_from_zenodo.sh` script
INPUT_PATH=`pwd`/../images/10.5281_zenodo.7057076
# Define a unique output path that depends on the username
OUTPUT_PATH=/data/active/jluethi/Fractal/20221004_2x2_Zenodo_3D
rm -rv $OUTPUT_PATH

TMPJSON=${TMPDIR}/tmp.json
TMPTASKS=${TMPDIR}/core_tasks.json

CMD="fractal"
CMD_JSON="python aux_extract_from_simple_json.py $TMPJSON"
fractal task list > $TMPTASKS

# Create project
$CMD -j project new $PRJ_NAME $TMPDIR > $TMPJSON
PRJ_ID=`$CMD_JSON project_id`
DS_IN_ID=`$CMD_JSON dataset_id "default"`
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

echo "{\"num_levels\": 5, \"executor\": \"cpu-low\", \"coarsening_xy\": 2, \"channel_parameters\": {\"A01_C01\": {\"label\": \"DAPI\",\"colormap\": \"00FFFF\",\"start\": 50,\"end\": 700 }, \"A01_C02\": {\"label\": \"nanog\",\"colormap\": \"FF00FF\",\"start\": 20,\"end\": 200 }, \"A02_C03\": {\"label\": \"Lamin B1\",\"colormap\": \"FFFF00\",\"start\": 50,\"end\": 1500 }}}" > ${TMPDIR}/args_create.json
fractal task add-subtask $WF_ID "Create OME-ZARR structure" --args-file ${TMPDIR}/args_create.json

echo "{\"executor\": \"cpu-low\"}" > ${TMPDIR}/args_yoko.json
fractal task add-subtask $WF_ID "Yokogawa to Zarr" > ${TMPDIR}/args_yoko.json

echo "{\"labeling_level\": 1, \"executor\": \"gpu\", \"ROI_table_name\": \"FOV_ROI_table\"}" > ${TMPDIR}/args_labeling_FOV.json
fractal task add-subtask $WF_ID "Per-FOV image labeling" --args-file ${TMPDIR}/args_labeling_FOV.json

echo "{\"level\": 0, \"measurement_table_name\": \"nuclei_3D\", \"executor\": \"cpu-high\", \"ROI_table_name\": \"FOV_ROI_table\",\"workflow_file\": \"$TMPDIR/../regionprops_from_existing_labels_feature.yaml\"}" > ${TMPDIR}/args_measurement_FOV.json
fractal task add-subtask $WF_ID "Measurement" --args-file ${TMPDIR}/args_measurement_FOV.json

# Paths of illumination correction images need to be accessible on the server.
# This works if one runs the client from the same machine as the server. Otherwise, change `root_path_corr`
echo "{\"overwrite\": true, \"executor\": \"cpu-low\", \"dict_corr\": {\"root_path_corr\": \"$TMPDIR/../../../illum_corr_images/\", \"A01_C01\": \"20220621_UZH_manual_illumcorr_40x_A01_C01.png\", \"A01_C02\": \"20220621_UZH_manual_illumcorr_40x_A01_C02.png\", \"A02_C03\": \"20220621_UZH_manual_illumcorr_40x_A02_C03.png\"}}" > ${TMPDIR}/args_illum.json
fractal task add-subtask $WF_ID "Illumination correction" --args-file ${TMPDIR}/args_illum.json

echo "{\"executor\": \"cpu-low\"}" > ${TMPDIR}/args_replicate.json
fractal task add-subtask $WF_ID "Replicate Zarr structure" --args-file ${TMPDIR}/args_replicate.json

echo "{\"executor\": \"cpu-low\"}" > ${TMPDIR}/args_mip.json
fractal task add-subtask $WF_ID "Maximum Intensity Projection" --args-file ${TMPDIR}/args_mip.json

echo "{\"labeling_level\": 2, \"executor\": \"cpu-low\", \"ROI_table_name\": \"well_ROI_table\"}" > ${TMPDIR}/args_labeling.json
fractal task add-subtask $WF_ID "Per-FOV image labeling" --args-file ${TMPDIR}/args_labeling.json

echo "{\"level\": 0, \"measurement_table_name\": \"nuclei\", \"executor\": \"cpu-low\", \"ROI_table_name\": \"well_ROI_table\",\"workflow_file\": \"$TMPDIR/../regionprops_from_existing_labels_feature.yaml\"}" > ${TMPDIR}/args_measurement.json
fractal task add-subtask $WF_ID "Measurement" --args-file ${TMPDIR}/args_measurement.json

# Apply workflow
fractal task apply $PRJ_ID $DS_IN_ID $DS_OUT_ID $WF_ID

