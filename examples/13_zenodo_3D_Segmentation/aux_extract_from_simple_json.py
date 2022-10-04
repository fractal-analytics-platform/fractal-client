import json
import sys


msg = (
    "Usage 1: python aux_extract_id_from_project_json.py FILENAME"
    "project_id\nUsage 2: python aux_extract_id_from_project_json.py"
    " FILENAME dataset_id DATASET_NAME"
)
if len(sys.argv[1:]) < 2:
    raise Exception(msg)

filename, which_id = sys.argv[1:3]
if which_id == "dataset_id":
    if len(sys.argv[1:]) != 3:
        raise Exception(msg)
    dataset_name = sys.argv[3]

with open(filename, "r") as fin:
    d = json.load(fin)

if which_id == "project_id":
    # Safety check
    print(d["id"])
elif which_id == "dataset_id":
    dataset_list = d["dataset_list"]
    dataset_id = next(
        ds["id"] for ds in dataset_list if ds["name"] == dataset_name
    )
    print(dataset_id)
else:
    raise Exception("ERROR: {which_id=}")