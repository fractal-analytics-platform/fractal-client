import os
import sys

from fractal_cmd import dataset_update_type
from fractal_cmd import datasets_add_resources
from fractal_cmd import datasets_list
from fractal_cmd import project_new
from fractal_cmd import projects_list
from fractal_cmd import task_add
from fractal_cmd import task_list
from fractal_cmd import workflow_apply
from fractal_cmd import workflow_list
from fractal_cmd import workflow_new


# General variables and paths, to be set by user
resource_in = (
    "/data/active/fractal-temp/3D/PelkmansLab/"
    + "CardiacMultiplexing/Cycle1_subset_extras/"
)
tmp_path = os.getcwd() + "/Temporary_folder_for_tests/"
resource_out = tmp_path + "result/"


####################################
# Do not modify code under this line

# Quick&dirty way to ignore function decorators
# (which are otherwise used for CLI)
project_new = project_new.__dict__["callback"]
projects_list = projects_list.__dict__["callback"]
datasets_add_resources = datasets_add_resources.__dict__["callback"]
dataset_update_type = dataset_update_type.__dict__["callback"]
datasets_list = datasets_list.__dict__["callback"]
task_add = task_add.__dict__["callback"]
task_list = task_list.__dict__["callback"]
workflow_new = workflow_new.__dict__["callback"]
workflow_list = workflow_list.__dict__["callback"]
workflow_apply = workflow_apply.__dict__["callback"]


# Check that tmp_path folder is not there
if os.path.isdir(tmp_path):
    sys.exit(
        f"ERROR: {tmp_path} folder already exists, "
        f"you should remove it or change the tmp_path variable."
    )

# General variables and path, part 2 (do not modify this)
project_name = "mwe-test"
dataset_name = "dstest"
workflow_name = "wftest"
ds_type = "tif"
task_name = "yokogawa_to_zarr"
input_type = "tif"
output_type = "zarr"

# Prepare and execute a workflow
project_new(project_name, tmp_path, dataset_name)
datasets_add_resources(project_name, dataset_name, [resource_in])
dataset_update_type(project_name, dataset_name, ds_type)
task_add(task_name, input_type, output_type)
workflow_new(project_name, workflow_name, [task_name])
workflow_list(project_name)
# workflow_apply(project_name, workflow_name, [dataset_name], [dataset_name],
#                [resource_in], [resource_out])
