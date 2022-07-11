import json
import pathlib
import shutil

import pytest

from fractal.tasks.illumination_correction import correct
from fractal.tasks.illumination_correction import illumination_correction


@pytest.mark.skip(reason="This test is not ready yet")
def test_illumination_correction_number_executions(
    tmp_path: pathlib.Path, monkeypatch
):
    """
    GIVEN a zarr array stored on disk
    WHEN I apply illumination_correction
    THEN correct(..) is executed as many times as the number of images
    """

    # Copy a reference zarr into a temporary folder
    zarrurl = tmp_path.resolve().as_posix() + "/plate.zarr/"
    newzarrurl = tmp_path.resolve().as_posix() + "/plate_new.zarr/"
    shutil.copytree("tests/data/plate_ones.zarr", zarrurl)
    component = "B/03/0/"
    zarrurl += component
    newzarrurl += component

    # Load some useful variables
    with open(zarrurl + "0/.zarray", "r") as f:
        zarray = json.load(f)
    num_C = zarray["shape"][0] // zarray["chunks"][0]
    num_Z = zarray["shape"][1] // zarray["chunks"][1]
    num_Y = zarray["shape"][2] // zarray["chunks"][2]
    num_X = zarray["shape"][3] // zarray["chunks"][3]
    num_chunks = num_C * num_Z * num_Y * num_X

    # Patch correct() function, to keep track of the number of calls
    logfile = tmp_path.resolve().as_posix() + "log_function_correct.txt"
    logfile = "log_function_correct.txt"
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
    illumination_correction(
        zarrurl,
        overwrite=False,
        newzarrurl=newzarrurl,
        chl_list=["A01_C01", "A02_C02"],
        path_dict_corr="tests/data/illumination_correction/illum.json",
        coarsening_xy=2,
        background=100,
    )

    # Verify the total number of calls
    with open(logfile, "r") as f:
        tot_calls_correct = len(f.read().splitlines())
    assert tot_calls_correct == num_chunks


# ITERATIONS
# overwrite yes (existing or not) vs false
# Verify output
# the right number of times (equal to the number of blocks), and output must
# be of the right size
