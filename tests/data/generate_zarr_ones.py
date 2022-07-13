import json

import dask.array as da
import numpy as np


num_C = 2
num_Z = 1
num_Y = 2
num_X = 1
num_levels = 5

x = da.ones(
    (num_C, num_Z, num_Y * 2160, num_X * 2560), chunks=(1, 1, 2160, 2560)
).astype(np.uint16)


zarrurl = "plate_ones.zarr/"
component = "B/03/0/"

for ind_level in range(num_levels):
    scale = 2**ind_level
    y = da.coarsen(np.min, x, {2: scale, 3: scale}).rechunk(
        (1, 1, 2160, 2560), balance=True
    )
    y.to_zarr(
        zarrurl, component=f"{component}{ind_level}", dimension_separator="/"
    )

zattrs = {
    "multiscales": [
        {
            "axes": [],
            "datasets": [{"path": level} for level in range(num_levels)],
            "version": "0.3",
        }
    ]
}
with open(f"{zarrurl}{component}.zattrs", "w") as jsonfile:
    json.dump(zattrs, jsonfile)
