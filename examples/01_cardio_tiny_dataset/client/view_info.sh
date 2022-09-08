# NOTE: This works if env variables are already defined,
#       and if the user is already registered

CMD="poetry run fractal"

$CMD project list
$CMD task list
$CMD dataset show $PROJECT_NAME $DATASET_IN_NAME
$CMD dataset show-resources $PROJECT_NAME $DATASET_IN_NAME
$CMD dataset show $PROJECT_NAME $DATASET_OUT_NAME
$CMD dataset show-resources $PROJECT_NAME $DATASET_OUT_NAME
