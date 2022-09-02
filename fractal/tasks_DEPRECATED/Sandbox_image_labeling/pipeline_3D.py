import os
import shutil
from pathlib import Path

from devtools import debug

from fractal.tasks.create_zarr_structure import create_zarr_structure
from fractal.tasks.image_labeling import image_labeling
from fractal.tasks.yokogawa_to_zarr import yokogawa_to_zarr


# Set images/zarr folders
img_path = Path(
    "/home/tommaso/Fractal/mwe_fractal/fractal/tasks/"
    + "Sandbox_image_labeling/images/*.png"
)
zarr_path = Path("./tmp_out/*.zarr")
zarr_mip_path = Path("./tmp_out_mip/*.zarr")
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
num_levels = 5
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
        rows=2,
        cols=2,
        metadata=metadata,
        component=component,
    )
debug(metadata)

# Labeling
for component in metadata["well"]:
    image_labeling(
        input_paths=[zarr_path],
        output_path=zarr_path,
        component=component,
        metadata=metadata,
        labeling_channel="A01_C01",
        labeling_level=4,
        num_threads=1,
    )
