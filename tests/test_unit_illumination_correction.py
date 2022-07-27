import json
import pathlib
import shutil

import anndata as ad
import dask.array as da
import numpy as np
import pytest
from pytest import MonkeyPatch

from fractal.tasks.illumination_correction import correct
from fractal.tasks.illumination_correction import illumination_correction
from fractal.tasks.lib_regions_of_interest import convert_ROI_table_to_indices
from fractal.tasks.lib_regions_of_interest import (
    split_3D_indices_into_z_layers,
)
from fractal.tasks.lib_zattrs_utils import extract_zyx_pixel_sizes


@pytest.mark.parametrize("overwrite", [True, False])
def test_illumination_correction(
    overwrite: bool,
    tmp_path: pathlib.Path,
    monkeypatch: MonkeyPatch,
):
    # GIVEN a zarr pyramid on disk, made of all ones
    # WHEN I apply illumination_correction
    # THEN correct(..) is executed as many times as
    #      (number of FOVs) x (number of channels)
    # AND the output array has ones at all pyramid levels

    # Copy a reference zarr into a temporary folder
    # FIXME: replace the try/except with a relative path
    zarrurl = (tmp_path / "plate.zarr").resolve().as_posix()
    try:
        shutil.copytree("tests/data/plate_ones.zarr", zarrurl)
    except FileNotFoundError as e:
        raise FileNotFoundError(
            f"{e}\nProbably you are not running from fractal root folder"
        )
    zarrurl += "/B/03/0/"

    # Set newzarrurl
    if overwrite:
        newzarrurl = zarrurl
    else:
        newzarrurl = zarrurl.replace("plate.zarr", "newplate.zarr")

    # Read FOV ROIs and create corresponding indices
    pixels = extract_zyx_pixel_sizes(zarrurl + ".zattrs", level=0)
    ROIs = ad.read_zarr(zarrurl + "tables/FOV_ROI_table/")
    list_indices = convert_ROI_table_to_indices(
        ROIs, level=0, full_res_pxl_sizes_zyx=pixels
    )
    list_indices = split_3D_indices_into_z_layers(list_indices)
    num_FOVs = len(list_indices)

    # Load some useful variables
    with open(zarrurl + "0/.zarray", "r") as f:
        zarray = json.load(f)
    shape = zarray["shape"]
    num_channels = shape[0]

    # Patch correct() function, to keep track of the number of calls
    logfile = (tmp_path / "log_function_correct.txt").resolve().as_posix()
    with open(logfile, "w") as log:
        log.write("")

    def patched_correct(*args, **kwargs):
        with open(logfile, "a") as log:
            log.write("1\n")
        return correct(*args, **kwargs)

    monkeypatch.setattr(
        "fractal.tasks.illumination_correction.correct", patched_correct
    )

    # Call illumination correction task, with patched correct()
    # FIXME: make tests/data paths relative?
    illumination_correction(
        zarrurl,
        overwrite=overwrite,
        newzarrurl=newzarrurl,
        chl_list=["A01_C01", "A02_C02"],
        path_dict_corr="tests/data/illumination_correction/illum.json",
        coarsening_xy=2,
        background=0,
    )

    # FIXME: make tests/data paths relative?
    if not overwrite:
        shutil.copy("tests/data/plate_ones.zarr/B/03/0/.zattrs", newzarrurl)

    # Verify the total number of calls
    with open(logfile, "r") as f:
        tot_calls_correct = len(f.read().splitlines())
    assert tot_calls_correct == (num_channels * num_FOVs)

    # Verify the output
    num_levels = 5
    for ind_level in range(num_levels):
        old = da.from_zarr(f"tests/data/plate_ones.zarr/B/03/0/{ind_level}")
        new = da.from_zarr(f"{newzarrurl}{ind_level}")
        assert old.shape == new.shape
        assert old.chunks == new.chunks
        assert new.compute()[0, 0, 0, 0] == 1
        assert np.allclose(old.compute(), new.compute())
