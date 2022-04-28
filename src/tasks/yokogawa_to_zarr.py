import os
import re
from glob import glob

import dask.array as da
import numpy as np
from dask import delayed
from skimage.io import imread


def sort_fun(s):

    """
    sort_fun takes a string (filename of a yokogawa images),
    extract site and z-index metadata and returns them as a list.

    :param s: filename
    :type s: str

    """

    site = re.findall(r"F(.*)L", s)[0]
    zind = re.findall(r"Z(.*)C", s)[0]
    return [site, zind]


def yokogawa_to_zarr(
    in_path, out_path, zarrurl, delete_in, rows, cols, ext, chl_list
):

    """
    Convert Yokogawa output (png, tif) to zarr file

    :param in_path: directory containing the input files
    :type in_path: str
    :param out_path: directory containing the output files
    :type out_path: str
    :param zarrurl: structure of the zarr folder
    :type zarrurl: str
    :param delete_in: delete input files, and folder if empty
    :type delete_in: bool
    :param rows: number of rows of the plate
    :type rows: int
    :param cols: number of columns of the plate
    :type cols: int
    :param ext: source images extension
    :type ext: str
    :param chl_list: list of the channels
    :type chl_list: list

    """

    r = zarrurl.split("/")[1]
    c = zarrurl.split("/")[2]

    lazy_imread = delayed(imread)
    fc_list = []
    fc1_list = []
    fc2_list = []
    fc3_list = []
    fc4_list = []

    print(chl_list)

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

        for r in range(int(rows)):
            cell = []

            for c in range(int(cols)):
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
    from argparse import ArgumentParser

    parser = ArgumentParser(prog="Yokogawa_to_zarr")

    parser.add_argument(
        "-i", "--in_path", help="directory containing the input files"
    )

    parser.add_argument(
        "-o", "--out_path", help="directory containing the output files"
    )

    parser.add_argument(
        "-z",
        "--zarrurl",
        help="structure of the zarr folder",
    )

    parser.add_argument(
        "-d",
        "--delete_in",
        help="Delete input files and folder",
    )

    parser.add_argument(
        "-r",
        "--rows",
        help="Number of rows of final image",
    )

    parser.add_argument(
        "-c",
        "--cols",
        help="Number of columns of final image",
    )

    parser.add_argument(
        "-e",
        "--ext",
        help="source images extension",
    )

    parser.add_argument(
        "-C",
        "--chl_list",
        nargs="+",
        help="list of channels ",
    )

    args = parser.parse_args()

    yokogawa_to_zarr(
        args.in_path,
        args.out_path,
        args.zarrurl,
        args.delete_in,
        args.rows,
        args.cols,
        args.ext,
        args.chl_list,
    )
