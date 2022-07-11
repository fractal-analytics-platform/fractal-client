import json

import dask.array as da
import numpy as np

num_C = 2
num_Z = 1
num_Y = 2
num_X = 1
x = da.ones(
    (num_C, num_Z, num_Y * 2160, num_X * 2560), chunks=(1, 1, 2160, 2560)
).astype(np.uint16)


zarrurl = "plate_ones.zarr/"
component = "B/03/0/"
x.to_zarr(zarrurl, component=f"{component}{0}", dimension_separator="/")


zattrs = {
    "multiscales": [
        {
            "axes": [],
            "datasets": [{"path": level} for level in range(5)],
            "version": "0.3",
        }
    ]
}
with open(f"{zarrurl}{component}.zattrs", "w") as jsonfile:
    json.dump(zattrs, jsonfile)
