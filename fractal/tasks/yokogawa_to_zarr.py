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
    zarrurl,
    in_path=None,
    ext=None,
    dims=None,
    chl_list=None,
    num_levels=5,
    coarsening_factor_xy=2,
    coarsening_factor_z=1,
    delete_in=False,
):

    """
    Convert Yokogawa output (png, tif) to zarr file

    #FIXME docstring

    :param in_path: directory containing the input files
    :type in_path: str
    :param zarrurl: structure of the zarr folder
    :type zarrurl: str
    :param delete_in: delete input files, and folder if empty
    :type delete_in: str
    :param rows: number of rows of the well
    :type rows: int
    :param cols: number of columns of the well
    :type cols: int
    :param ext: source images extension
    :type ext: str
    :param chl_list: list of the channels
    :type chl_list: list
    :param num_levels: number of levels in the zarr pyramid
    :type num_levels: int
    :param coarsening_factor_xy: coarsening factor along X,Y
    :type coarsening_factor_xy: int
    :param coarsening_factor_z: coarsening factor along Z
    :type coarsening_factor_z: int

    """

    # Hard-coded values (by now) of how chunk size to be passed to rechunk,
    # both at level 0 (before coarsening) and at levels 1,2,.. (after
    # repeated coarsening). Note that balance=True may override these values.
    chunk_size_x = 256 * 5
    chunk_size_y = 216 * 5

    # Define well
    if not zarrurl.endswith("/"):
        zarrurl += "/"
    r = zarrurl.split("/")[-4]
    c = zarrurl.split("/")[-3]

    # Define grid of sites (within the well)
    rows, cols = dims[:]

    lazy_imread = delayed(imread)
    fc_list = {level: [] for level in range(num_levels)}

    print(chl_list)

    for ch in chl_list:

        l_rows = []
        all_rows = []

        print(zarrurl, r, c)
        filenames = sorted(
            glob(in_path + f"*_{r+c}_*C{ch}." + ext), key=sort_fun
        )
        print(in_path + f"*_{r+c}_*C{ch}." + ext)
        max_z = max([re.findall(r"Z(.*)C", s)[0] for s in filenames])

        sample = imread(filenames[0])

        s = 0
        e = int(max_z)

        for ro in range(int(rows)):
            cell = []

            for co in range(int(cols)):
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
        # At this point, all_rows has four dimensions: z,site,y,x

        # Define coarsening options
        coarsening = {1: coarsening_factor_xy, 2: coarsening_factor_xy}
        f_matrices = {}
        for level in range(num_levels):
            if level == 0:
                f_matrices[level] = da.concatenate(all_rows, axis=1).rechunk(
                    {1: chunk_size_x, 2: chunk_size_y}, balance=True
                )
                # After concatenate, f_matrices[0] has three dimensions: z,y,x
                if coarsening_factor_z > 1:
                    f_matrices[level] = da.coarsen(
                        np.min,
                        f_matrices[level],
                        {0: coarsening_factor_z},
                        trim_excess=True,
                    )
            else:
                # FIXME: add here a threshold to avoid having a lot
                # of minuscule chunks
                f_matrices[level] = da.coarsen(
                    np.min,
                    f_matrices[level - 1],
                    coarsening,
                    trim_excess=True,
                ).rechunk(
                    {
                        1: max(
                            1, chunk_size_x // (coarsening_factor_xy**level)
                        ),
                        2: max(
                            1, chunk_size_y // (coarsening_factor_xy**level)
                        ),
                    },
                    balance=True,
                )
            fc_list[level].append(f_matrices[level])

    level_data = [
        da.stack(fc_list[level], axis=0) for level in range(num_levels)
    ]

    shape_list = []
    for i, level in enumerate(level_data):
        level.to_zarr(zarrurl + f"{i}/", dimension_separator="/")
        print(f"Chunks at level {i}:\n", level.chunks)
        shape_list.append(level.shape)
    print()

    if delete_in == "True":
        for f in filenames:
            try:
                os.remove(f)
            except OSError as e:
                print("Error: %s : %s" % (f, e.strerror))
    return shape_list


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(prog="Yokogawa_to_zarr")

    parser.add_argument(
        "-i", "--in_path", help="directory containing the input files"
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

    parser.add_argument(
        "-nl",
        "--num_levels",
        type=int,
        help="number of levels in the Zarr pyramid",
    )

    parser.add_argument(
        "-cxy",
        "--coarsening_xy",
        default=2,
        type=int,
        help="coarsening factor along X and Y (optional, defaults to 2)",
    )

    parser.add_argument(
        "-cz",
        "--coarsening_z",
        default=1,
        type=int,
        help="coarsening factor along Z (optional, defaults to 1)",
    )

    args = parser.parse_args()

    yokogawa_to_zarr(
        args.in_path,
        args.zarrurl,
        args.delete_in,
        args.rows,
        args.cols,
        args.ext,
        args.chl_list,
        args.num_levels,
        args.coarsening_xy,
        args.coarsening_z,
    )
