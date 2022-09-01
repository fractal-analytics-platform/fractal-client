# NOTE: This works if env variables are already defined,
#       and if the user is already registered

poetry run client project list
poetry run client task list
poetry run client dataset show $PROJECT_NAME $DATASET_IN_NAME
poetry run client dataset show-resources $PROJECT_NAME $DATASET_IN_NAME
poetry run client dataset show $PROJECT_NAME $DATASET_OUT_NAME
poetry run client dataset show-resources $PROJECT_NAME $DATASET_OUT_NAME
