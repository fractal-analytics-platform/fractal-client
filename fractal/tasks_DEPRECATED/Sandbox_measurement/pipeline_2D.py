import os
import shutil
from pathlib import Path

from devtools import debug

from fractal.tasks.create_zarr_structure import create_zarr_structure
from fractal.tasks.image_labeling_whole_well import image_labeling_whole_well
from fractal.tasks.maximum_intensity_projection import (
    maximum_intensity_projection,
)
from fractal.tasks.measurement import measurement
from fractal.tasks.replicate_zarr_structure import replicate_zarr_structure
from fractal.tasks.yokogawa_to_zarr import yokogawa_to_zarr


# Set images/zarr folders
img_path = Path("images/*.png")
zarr_path = Path("tmp_out/*.zarr")
zarr_mip_path = Path("tmp_out_mip/*.zarr")
if os.path.isdir(zarr_path.parent):
    shutil.rmtree(zarr_path.parent)
if os.path.isdir(zarr_mip_path.parent):
    shutil.rmtree(zarr_mip_path.parent)

# Set useful parameters
channel_parameters = {
    "A01_C01": {
        "label": "DAPI",
        "colormap": "00FFFF",
        "start": 0,
        "end": 700,
    },
    "A01_C02": {
        "label": "nanog",
        "colormap": "FF00FF",
        "start": 0,
        "end": 180,
    },
    "A02_C03": {
        "label": "Lamin B1",
        "colormap": "FFFF00",
        "start": 0,
        "end": 1500,
    },
}
num_levels = 4
coarsening_xy = 2


# START WORFKLOW

metadata = {}

# Create zarr structure
metadata_update = create_zarr_structure(
    input_paths=[img_path],
    output_path=zarr_path,
    channel_parameters=channel_parameters,
    num_levels=num_levels,
    coarsening_xy=coarsening_xy,
    metadata_table="mrf_mlf",
)
metadata.update(metadata_update)
debug(metadata)

# Yokogawa to zarr
for component in metadata["well"]:
    yokogawa_to_zarr(
        input_paths=[zarr_path],
        output_path=zarr_path,
        rows=1,
        cols=2,
        metadata=metadata,
        component=component,
    )
debug(metadata)

# Replicate zarr structure
metadata_update = replicate_zarr_structure(
    input_paths=[zarr_path],
    output_path=zarr_mip_path,
    suffix="mip",
    metadata=metadata,
)
metadata.update(metadata_update)
debug(metadata)

# Maximum intensity projection
for component in metadata["well"]:
    maximum_intensity_projection(
        input_paths=[zarr_mip_path],
        output_path=zarr_mip_path,
        component=component,
        metadata=metadata,
    )


# Labeling
for component in metadata["well"]:
    image_labeling_whole_well(
        input_paths=[zarr_mip_path],
        output_path=zarr_mip_path,
        component=component,
        metadata=metadata,
        labeling_channel="A01_C01",
        labeling_level=2,
    )

# Measurement
for component in metadata["well"]:
    measurement(
        input_paths=[zarr_mip_path],
        output_path=zarr_mip_path,
        metadata=metadata,
        component=component,
        labeling_channel="A01_C01",
        level=0,
        workflow_file="regionprops_from_existing_labels.yaml",
        table_name="nuclei",
        whole_well=True,
    )
