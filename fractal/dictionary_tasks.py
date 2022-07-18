"""
Copyright 2022 (C) Friedrich Miescher Institute for Biomedical Research and
University of Zurich

Original authors:
Tommaso Comparin <tommaso.comparin@exact-lab.it>
Marco Franzon <marco.franzon@exact-lab.it>

This file is part of Fractal and was originally developed by eXact lab S.r.l.
<exact-lab.it> under contract with Liberali Lab from the Friedrich Miescher
Institute for Biomedical Research and Pelkmans Lab from the University of
Zurich.
"""
from fractal.tasks.create_zarr_structure import create_zarr_structure
from fractal.tasks.create_zarr_structure_multifov import (
    create_zarr_structure_multifov,
)
from fractal.tasks.illumination_correction import illumination_correction
from fractal.tasks.image_labeling import image_labeling
from fractal.tasks.image_labeling_whole_well import image_labeling_whole_well
from fractal.tasks.maximum_intensity_projection import (
    maximum_intensity_projection,
)
from fractal.tasks.replicate_zarr_structure import replicate_zarr_structure
from fractal.tasks.replicate_zarr_structure_mip import (
    replicate_zarr_structure_mip,
)
from fractal.tasks.yokogawa_to_zarr import yokogawa_to_zarr
from fractal.tasks.yokogawa_to_zarr_multifov import yokogawa_to_zarr_multifov

dict_tasks = {}
dict_tasks["create_zarr_structure"] = create_zarr_structure
dict_tasks["create_zarr_structure_multifov"] = create_zarr_structure_multifov
dict_tasks["yokogawa_to_zarr"] = yokogawa_to_zarr
dict_tasks["yokogawa_to_zarr_multifov"] = yokogawa_to_zarr_multifov
dict_tasks["replicate_zarr_structure_mip"] = replicate_zarr_structure_mip
dict_tasks["maximum_intensity_projection"] = maximum_intensity_projection
dict_tasks["replicate_zarr_structure"] = replicate_zarr_structure
dict_tasks["illumination_correction"] = illumination_correction
dict_tasks["image_labeling"] = image_labeling
dict_tasks["image_labeling_whole_well"] = image_labeling_whole_well
