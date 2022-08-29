import json
import os

import anndata as ad
import dask.array as da
import numpy as np
import pandas as pd
import zarr
from anndata.experimental import write_elem
from napari_workflows._io_yaml_v1 import load_workflow

from fractal.tasks.lib_regions_of_interest import convert_ROI_table_to_indices
from fractal.tasks.lib_zattrs_utils import extract_zyx_pixel_sizes


# Inputs
in_path = "/data/active/fractal/tests/Temporary_data_UZH_1_well_2x2_sites"
component = "20200812-CardiomyocyteDifferentiation14-Cycle1.zarr/B/03/0"
labeling_channel = "A01_C01"
level = 0
metadata = dict(coarsening_xy=2, channel_list=["A01_C01"])
workflow_file = "regionprops_from_existing_labels.yaml"
table_name = "nuclei"


# Load metadata
coarsening_xy = metadata["coarsening_xy"]
chl_list = metadata["channel_list"]

# Find channel index
if labeling_channel not in chl_list:
    raise Exception(f"ERROR: {labeling_channel} not in {chl_list}")
ind_channel = chl_list.index(labeling_channel)

# Load zattrs file
zattrs_file = f"{in_path}/{component}/.zattrs"
with open(zattrs_file, "r") as jsonfile:
    zattrs = json.load(jsonfile)

# Try to read channel label from OMERO metadata
try:
    omero_label = zattrs["omero"]["channels"][ind_channel]["label"]
    label_name = f"label_{omero_label}"
except (KeyError, IndexError):
    label_name = f"label_{ind_channel}"

# Set level=0, to avoid possible errors, see
# https://github.com/fractal-analytics-platform/fractal/issues/69#issuecomment-1230074703
if level > 0:
    raise Exception("Measurement should be for level=0")

# Load
img = da.from_zarr(f"{in_path}/{component}/{level}")[ind_channel]
label_img = da.from_zarr(f"{in_path}/{component}/labels/{label_name}/{level}")

# Find upscale_factor for labels array
img_shape = img.shape
label_shape = label_img.shape
upscale_factor = img_shape[1] // label_shape[1]
if img_shape[0] != label_shape[0]:
    raise Exception("Error in dapi/label array shapes", img_shape, label_shape)

if (
    label_shape[1] * upscale_factor != img_shape[1]
    or upscale_factor != img_shape[2] // label_shape[2]
):
    raise Exception(
        "Error in dapi/label array shapes",
        img_shape,
        label_shape,
        f"with {upscale_factor=}",
    )

# Upscale labels array - see https://stackoverflow.com/a/7525345/19085332
label_img_up_x = np.repeat(label_img, upscale_factor, axis=2)
label_img_up = np.repeat(label_img_up_x, upscale_factor, axis=1)
if not label_img_up.shape == img_shape:
    raise Exception(
        "Error in dapi/label array shapes (after upscaling)",
        img.shape,
        label_img_up.shape,
    )

# Read FOV ROIs
FOV_ROI_table = ad.read_zarr(f"{in_path}/{component}/tables/FOV_ROI_table")

# Read pixel sizes from zattrs file
full_res_pxl_sizes_zyx = extract_zyx_pixel_sizes(
    f"{in_path}/{component}/.zattrs", level=0
)

# Create list of indices for 3D FOVs spanning the entire Z direction
list_indices = convert_ROI_table_to_indices(
    FOV_ROI_table,
    level=level,
    coarsening_xy=coarsening_xy,
    full_res_pxl_sizes_zyx=full_res_pxl_sizes_zyx,
)

# Check that the target group is not already there, or fail fast
target_table_folder = f"{in_path}/{component}/tables/{table_name}"
if os.path.isdir(target_table_folder):
    raise Exception(f"ERROR: {target_table_folder} already exists.")


# Get the workflow
napari_workflow = load_workflow(workflow_file)

print(f"Workflow file:         {workflow_file}")
print(f"Resolution level:      {level}")
print(f"Labels upscale factor: {upscale_factor}")
print(f"Whole-array shape:     {img_shape}")

# Loop over FOV ROIs
list_dfs = []
for indices in list_indices:
    s_z, e_z, s_y, e_y, s_x, e_x = indices[:]
    ROI = (slice(s_z, e_z), slice(s_y, e_y), slice(s_x, e_x))
    print(f"Single-ROI shape:      {img[ROI].shape}")

    # Set the input images: DAPI channel image & label image for current ROI
    napari_workflow.set("dapi_img", img[ROI])
    napari_workflow.set("label_img", label_img_up[ROI])

    # Run the workflow
    df = napari_workflow.get("regionprops_DAPI")
    list_dfs.append(df)

# Concatenate all FOV dataframes
df_well = pd.concat(list_dfs, axis=0)

# Convert all to float (warning: some would be int, in principle)
measurement_dtype = np.float32
df_well = df_well.astype(measurement_dtype)

# Drop unnecessary columns
df_well.drop(labels=["label"], axis=1, inplace=True)

# Reset index (to avoid non-unique labels), and set it to str
# (to avoid ImplicitModificationWarning in anndata)
df_well.index = np.arange(0, len(df_well.index)).astype(str)

# Convert to anndata
measurement_table = ad.AnnData(df_well, dtype=measurement_dtype)

# Write to zarr group
group_tables = zarr.group(f"{in_path}/{component}/tables/")
write_elem(group_tables, table_name, measurement_table)
