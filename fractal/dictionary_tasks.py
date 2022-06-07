from fractal.tasks.create_zarr_structure import create_zarr_structure
from fractal.tasks.illumination_correction import illumination_correction
from fractal.tasks.maximum_intensity_projection import (
    maximum_intensity_projection,
)
from fractal.tasks.replicate_zarr_structure import replicate_zarr_structure
from fractal.tasks.replicate_zarr_structure_mip import (
    replicate_zarr_structure_mip,
)
from fractal.tasks.yokogawa_to_zarr import yokogawa_to_zarr

dict_tasks = {}
dict_tasks["create_zarr_structure"] = create_zarr_structure
dict_tasks["yokogawa_to_zarr"] = yokogawa_to_zarr
dict_tasks["replicate_zarr_structure_mip"] = replicate_zarr_structure_mip
dict_tasks["maximum_intensity_projection"] = maximum_intensity_projection
dict_tasks["replicate_zarr_structure"] = replicate_zarr_structure
dict_tasks["illumination_correction"] = illumination_correction
