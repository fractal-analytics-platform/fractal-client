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
import numpy as np
import pytest

from fractal.tasks.yokogawa_to_zarr import yokogawa_to_zarr


f1 = "plate_well_T0001F001L01A01Z01C01.png"
f2 = "plate_well_T0001F002L01A01Z01C01.png"
f3 = "plate_well_T0001F003L01A01Z01C01.png"
f4 = "plate_well_T0001F004L01A01Z01C01.png"

chl_list = ["01"]
num_levels = 5
coarsening_factor_xy = 2
coarsening_factor_z = 1


@pytest.mark.skip(
    reason="TODO: update after porting to new server-based architecture"
)
def test_yokogawa_to_zarr(mocker):

    mocker.patch(
        "fractal.tasks.yokogawa_to_zarr.sorted", return_value=[f1, f2, f3, f4]
    )

    mocker.patch("fractal.tasks.yokogawa_to_zarr.max", return_value="01")

    mocker.patch(
        "fractal.tasks.yokogawa_to_zarr.imread", return_value=np.ones((64, 64))
    )

    mocker.patch(
        "dask.delayed",
        new_callable=mocker.PropertyMock,
        return_value=[
            np.ones((64, 64)),
            np.ones((64, 64)),
            np.ones((64, 64)),
            np.ones((64, 64)),
        ],
    )

    mocker.patch(
        "dask.array.core.to_zarr", return_value=np.ones((64, 64)).shape
    )

    res = yokogawa_to_zarr(
        "plate.zarr/row/column/fov/",
        in_path="/tmp/",
        ext="png",
        rows=2,
        cols=2,
        chl_list=["A01_C01"],
        num_levels=5,
        coarsening_xy=2,
    )

    assert res == [
        (1, 1, 128, 128),
        (1, 1, 64, 64),
        (1, 1, 32, 32),
        (1, 1, 16, 16),
        (1, 1, 8, 8),
    ]
