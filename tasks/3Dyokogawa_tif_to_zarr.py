import os
import re
import sys
from glob import glob

import dask.array as da
import numpy as np
from dask import delayed
from skimage.io import imread


def sort_fun(s):
    site = re.findall(r"F(.*)L", s)[0]
    zind = re.findall(r"Z(.*)C", s)[0]
    return [site, zind]


def yokogawa_tif_to_zarr(
    in_path, out_path, zarrurl, delete_in, rows, cols, ext, chl_list
):

    r = zarrurl.split("/")[1]
    c = zarrurl.split("/")[2]

    lazy_imread = delayed(imread)
    fc_list = []
    fc1_list = []
    fc2_list = []
    fc3_list = []
    fc4_list = []

    for ch in chl_list:

        l_rows = []
        all_rows = []

        filenames = sorted(
            glob(in_path + f"*_{r+c}_*C{ch}*." + ext), key=sort_fun
        )
        print(in_path + f"*_{r+c}_*C{ch}*." + ext)
        max_z = max([re.findall(r"Z(.*)C", s)[0] for s in filenames])

        sample = imread(filenames[0])

        s = 0
        e = int(max_z)

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
        f1_matrix = [da.coarsen(np.min, x, {0: 2, 1: 2}) for x in f_matrix]
        f2_matrix = [da.coarsen(np.min, x, {0: 2, 1: 2}) for x in f1_matrix]
        f3_matrix = [da.coarsen(np.min, x, {0: 2, 1: 2}) for x in f2_matrix]
        f4_matrix = [da.coarsen(np.min, x, {0: 2, 1: 2}) for x in f3_matrix]

        fc_list.append(f_matrix)
        fc1_list.append(f1_matrix)
        fc2_list.append(f2_matrix)
        fc3_list.append(f3_matrix)
        fc4_list.append(f4_matrix)

    fc_stack = da.stack(fc_list, axis=0)
    fc1_stack = da.stack(fc1_list, axis=0)
    fc2_stack = da.stack(fc2_list, axis=0)
    fc3_stack = da.stack(fc3_list, axis=0)
    fc4_stack = da.stack(fc4_list, axis=0)

    tmp_lvl = [fc_stack, fc1_stack, fc2_stack, fc3_stack, fc4_stack]

    for i, level in enumerate(tmp_lvl):
        level.to_zarr(out_path + zarrurl + f"{i}/", dimension_separator="/")

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
    chl_list = sys.argv[8:]

    yokogawa_tif_to_zarr(
        in_path, out_path, zarrurl, delete_in, rows, cols, ext, chl_list
    )
