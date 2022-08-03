"""
Copyright 2022 (C) Friedrich Miescher Institute for Biomedical Research and
University of Zurich

Original authors:
Marco Franzon <marco.franzon@exact-lab.it>
Tommaso Comparin <tommaso.comparin@exact-lab.it>

This file is part of Fractal and was originally developed by eXact lab S.r.l.
<exact-lab.it> under contract with Liberali Lab from the Friedrich Miescher
Institute for Biomedical Research and Pelkmans Lab from the University of
Zurich.
"""
import pathlib

import numpy as np
from devtools import debug
from pytest import MonkeyPatch

from fractal.tasks.yokogawa_to_zarr import yokogawa_to_zarr


images = [
    "plate_well_T0001F001L01A01Z01C01.png",
    "plate_well_T0001F002L01A01Z01C01.png",
    "plate_well_T0001F003L01A01Z01C01.png",
    "plate_well_T0001F004L01A01Z01C01.png",
]

chl_list = ["01"]
num_levels = 5
coarsening_factor_xy = 2
coarsening_factor_z = 1


def test_yokogawa_to_zarr(
    mocker,
    tmp_path: pathlib.Path,
    monkeypatch: MonkeyPatch,
):

    debug(tmp_path)

    # Mock list of images
    mocker.patch("fractal.tasks.yokogawa_to_zarr.sorted", return_value=images)

    # Mock maximum Z-plane index
    mocker.patch("fractal.tasks.yokogawa_to_zarr.max", return_value="01")

    # Patch correct() function, to keep track of the number of calls
    logfile = (tmp_path / "log_function_correct.txt").resolve().as_posix()
    with open(logfile, "w") as log:
        log.write("")

    logfile = (tmp_path / "log_function_correct.txt").resolve().as_posix()

    def patched_imread(*args, **kwargs):
        with open(logfile, "a") as log:
            log.write("1\n")
        return np.ones((64, 64), dtype=np.uint16)

    monkeypatch.setattr(
        "fractal.tasks.yokogawa_to_zarr.imread", patched_imread
    )

    yokogawa_to_zarr(
        (tmp_path / "plate.zarr/row/column/fov/").as_posix(),
        in_path=tmp_path.as_posix(),
        ext="png",
        rows=2,
        cols=2,
        chl_list=["A01_C01"],
        num_levels=5,
        coarsening_xy=2,
    )

    # Read number of calls to imread
    num_calls_imread = np.loadtxt(logfile, dtype=int).sum()
    # Subtract one for each channel, for the dummy call at the beginning of
    # the task (used to determine shape and dtype)
    num_calls_imread -= len(chl_list)

    assert num_calls_imread == len(images)
