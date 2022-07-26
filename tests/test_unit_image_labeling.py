import json
import math
import pathlib
import shutil

import numpy as np
from pytest import MonkeyPatch

from fractal.tasks.image_labeling import image_labeling


def test_illumination_correction(
    tmp_path: pathlib.Path,
    monkeypatch: MonkeyPatch,
):
    # GIVEN a zarr pyramid on disk, made of all ones
    # WHEN I apply image_labeling
    # THEN segment_FOV(..) is executed as many times as the number of columns

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

    # Load some useful variables
    with open(zarrurl + "0/.zarray", "r") as f:
        zarray = json.load(f)
    shape = zarray["shape"]
    chunks = zarray["chunks"]
    num_columns = math.prod([shape[dim] // chunks[dim] for dim in [2, 3]])

    # Patch segment_FOV function, to keep track of the number of calls
    logfile = (tmp_path / "log_function_segment_FOV.txt").resolve().as_posix()
    with open(logfile, "w") as log:
        log.write("")

    def patched_segment_FOV(column, label_dtype=None, **kwargs):
        with open(logfile, "a") as log:
            log.write("1\n")
        return np.empty_like(column).astype(label_dtype)

    monkeypatch.setattr(
        "fractal.tasks.image_labeling.segment_FOV", patched_segment_FOV
    )

    image_labeling(
        zarrurl,
        coarsening_xy=2,
        labeling_level=0,
        labeling_channel="A01_C01",
        chl_list=["A01_C01", "A02_C02"],
        num_threads=2,
    )

    # Verify the total number of calls
    with open(logfile, "r") as f:
        tot_calls_correct = len(f.read().splitlines())
    assert tot_calls_correct == num_columns
