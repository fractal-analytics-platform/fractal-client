import numpy as np

from fractal.tasks.yokogawa_to_zarr import yokogawa_to_zarr


f1 = "plate_well_T0001F001L01A01Z01C01.png"
f2 = "plate_well_T0001F002L01A01Z01C01.png"
f3 = "plate_well_T0001F003L01A01Z01C01.png"
f4 = "plate_well_T0001F004L01A01Z01C01.png"

in_path = ""
out_path = ""
zarrurl = "plate.zarr/row/column/fov/"
delete_in = "False"
rows = "2"
cols = "2"
ext = "png"
chl_list = ["01"]
num_levels = 5
coarsening_factor_xy = 2
coarsening_factor_z = 1


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
        in_path,
        out_path,
        zarrurl,
        delete_in,
        rows,
        cols,
        ext,
        chl_list,
        num_levels,
        coarsening_factor_xy,
        coarsening_factor_z,
    )

    assert res == [
        (1, 1, 128, 128),
        (1, 1, 64, 64),
        (1, 1, 32, 32),
        (1, 1, 16, 16),
        (1, 1, 8, 8),
    ]
