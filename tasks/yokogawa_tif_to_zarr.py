import os
import sys
from glob import glob

import dask.array as da
from dask import delayed
from skimage.io import imread


# TODO use kwargs not args
def yokogawa_tif_to_zarr(
    in_path, zarrurl, zarr_f, delete_in, z_ind, cols, rows, ext
):

    plate = zarr_f.split("/")[0][:-5]
    well = zarr_f.split("/")[1]
    tims = zarr_f.split("/")[2]
    channel = zarr_f.split("/")[3]

    filenames = sorted(
        glob(in_path + f"{plate}_{well}_T{tims}F*Z{z_ind}C{channel}*." + ext)
    )
    # assuming that all files have same shape and type
    sample = imread(filenames[0])

    lazy_imread = delayed(imread)  # lazy reader
    lazy_arrays = [lazy_imread(fn) for fn in filenames]
    dask_arrays = [
        da.from_delayed(delayed_reader, shape=sample.shape, dtype=sample.dtype)
        for delayed_reader in lazy_arrays
    ]

    interval = len(dask_arrays) // int(rows)
    all_rows = []

    for j in range(int(rows)):
        t = j * interval
        tmp_row = []
        for i in range(int(cols)):
            tmp_row.append(dask_arrays[i + t])
        l_rows = da.concatenate(tmp_row, axis=1)
        all_rows.append(l_rows)
    f_matrix = da.concatenate(all_rows, axis=0)

    f_matrix.to_zarr(zarrurl + zarr_f + f"{z_ind}")

    if delete_in == "True":
        for f in filenames:
            try:
                os.remove(f)
            except OSError as e:
                print("Error: %s : %s" % (f, e.strerror))


if __name__ == "__main__":
    in_path = sys.argv[1]
    out_path = sys.argv[2]
    zarr_f = sys.argv[3]
    delete_in = sys.argv[4]
    cols = sys.argv[5]
    rows = sys.argv[6]
    ext = sys.argv[7]
    z_ind = sys.argv[8]

    yokogawa_tif_to_zarr(
        in_path, out_path, zarr_f, delete_in, z_ind, cols, rows, ext
    )
