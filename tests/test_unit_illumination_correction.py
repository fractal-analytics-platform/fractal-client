import json
import math
import pathlib
import shutil

import dask.array as da
import numpy as np
import pytest
from pytest import MonkeyPatch

from fractal.tasks.illumination_correction import correct
from fractal.tasks.illumination_correction import illumination_correction


@pytest.mark.parametrize("overwrite", [True, False])
@pytest.mark.skip("Missing ROIs table in test data")
def test_illumination_correction(
    overwrite: bool,
    tmp_path: pathlib.Path,
    monkeypatch: MonkeyPatch,
):
    # GIVEN a zarr pyramid on disk, made of all ones
    # WHEN I apply illumination_correction
    # THEN correct(..) is executed as many times as the number of chunks
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

    # Load some useful variables
    with open(zarrurl + "0/.zarray", "r") as f:
        zarray = json.load(f)
    shape = zarray["shape"]
    chunks = zarray["chunks"]
    num_chunks = math.prod([shape[dim] // chunks[dim] for dim in range(4)])

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
    assert tot_calls_correct == num_chunks

    # Verify the output
    num_levels = 5
    for ind_level in range(num_levels):
        old = da.from_zarr(f"tests/data/plate_ones.zarr/B/03/0/{ind_level}")
        new = da.from_zarr(f"{newzarrurl}{ind_level}")
        assert old.shape == new.shape
        assert old.chunks == new.chunks
        assert new.compute()[0, 0, 0, 0] == 1
        assert np.allclose(old.compute(), new.compute())
