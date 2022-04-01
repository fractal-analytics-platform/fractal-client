import os
import re
import sys
from glob import glob

import dask.array as da
from dask import delayed
from skimage.io import imread


def sort_fun(s):
    site = re.findall(r"F(.*)L", s)[0]
    zind = re.findall(r"Z(.*)C", s)[0]
    return [site, zind]


def yokogawa_tif_to_zarr(
    in_path, out_path, zarrurl, delete_in, rows, cols, ext, chl
):

    r = zarrurl.split("/")[1]
    c = zarrurl.split("/")[2]
    filenames = sorted(
        glob(in_path + f"*_{r+c}_*C{chl}*." + ext), key=sort_fun
    )

    max_z = max([re.findall(r"Z(.*)C", s)[0] for s in filenames])

    sample = imread(filenames[0])

    lazy_imread = delayed(imread)

    s = int(max_z)
    e = 0
    l_rows = []
    all_rows = []

    for j in range(int(rows)):
        cell = []
        for i in range(int(cols)):
            lazy_arrays = [lazy_imread(fn) for fn in filenames[s:e]]
            s += int(max_z)
            e += int(max_z)
            dask_arrays = [
                da.from_delayed(
                    delayed_reader, shape=sample.shape, dtype=sample.dtype
                )
                for delayed_reader in lazy_arrays
            ]
            z_stack = da.stack(dask_arrays, axis=0)
            cell.append(z_stack)
        l_rows = da.block(cell)
        all_rows.append(l_rows)
    f_matrix = da.concatenate(all_rows, axis=1)

    f_matrix.to_zarr(out_path + zarrurl, dimension_separator="/")

    if delete_in == "True":
        for f in filenames:
            try:
                os.remove(f)
            except OSError as e:
                print("Error: %s : %s" % (f, e.strerror))


if __name__ == "__main__":
    in_path = sys.argv[1]
    out_path = sys.argv[2]
    zarrurl = sys.argv[3]
    delete_in = sys.argv[4]
    rows = sys.argv[5]
    cols = sys.argv[6]
    ext = sys.argv[7]
    chl = sys.argv[8]

    yokogawa_tif_to_zarr(
        in_path, out_path, zarrurl, delete_in, rows, cols, ext, chl
    )
