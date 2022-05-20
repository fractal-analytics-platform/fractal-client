import os
import re
from glob import glob

import dask.array as da
import numpy as np
from dask import delayed
from skimage.io import imread


def sort_by_z(s):

    """
    sort_by_z takes a string (filename of a yokogawa images),
    extract z-index metadata and returns it.

    :param s: filename
    :type s: str

    """

    zind = re.findall(r"Z(.*)C", s)[0]
    return zind


def yokogawa_to_zarr(
    in_path,
    out_path,
    zarrurl,
    delete_in,
    sites_list,
    ext,
    chl_list,
    num_levels,
    coarsening_factor_xy,
    coarsening_factor_z,
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
    :type delete_in: str
    :param sites_list: list of sites
    :type sites_list: list
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

    r = zarrurl.split("/")[1]
    c = zarrurl.split("/")[2]

    lazy_imread = delayed(imread)

    print("channels:", chl_list)
    print("sites:", sites_list)

    # Loop over sites and channels
    for index_site, site in enumerate(sites_list):
        zarrurl_site = f"{zarrurl}{index_site}/"

        matrix_level_channel = {}

        for channel in chl_list:

            # Find all images for a given (site,channel)
            filenames = sorted(
                glob(in_path + f"*_{r+c}_*F{site}*C{channel}." + ext),
                key=sort_by_z,
            )
            print(in_path + f"*_{r+c}_*F{site}*C{channel}." + ext)

            # Build three-dimensional array for a given (site,channel)
            sample = imread(filenames[0])
            lazy_arrays = [lazy_imread(fn) for fn in filenames]
            dask_arrays = [
                da.from_delayed(
                    delayed_reader, shape=sample.shape, dtype=sample.dtype
                )
                for delayed_reader in lazy_arrays
            ]
            matrix_site_channel = da.stack(dask_arrays, axis=0)

            # Create pyramid
            coarsening = {
                1: coarsening_factor_xy,  # Y dimension
                2: coarsening_factor_xy,  # X dimension
            }
            for level in range(num_levels):
                if level == 0:
                    matrices_site_channel_levels = [matrix_site_channel]
                    if coarsening_factor_z > 1:
                        matrices_site_channel_levels[level] = da.coarsen(
                            np.min,
                            matrices_site_channel_levels[level],
                            {0: coarsening_factor_z},
                            trim_excess=True,
                        )
                else:
                    matrices_site_channel_levels.append(
                        da.coarsen(
                            np.min,
                            matrices_site_channel_levels[level - 1],
                            coarsening,
                            trim_excess=True,
                        )
                    )

                matrix_level_channel[
                    (level, channel)
                ] = matrices_site_channel_levels[level]

        for level in range(num_levels):
            level_list = [
                matrix_level_channel[key]
                for key in matrix_level_channel.keys()
                if key[0] == level
            ]
            level_stack = da.stack(level_list, axis=0)
            level_stack.to_zarr(
                out_path + zarrurl_site + f"{level}/", dimension_separator="/"
            )
            # shape_list.append(level.shape)

    shape_list = []

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
        "-S",
        "--sites_list",
        nargs="+",
        help="list of sites ",
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
        args.out_path,
        args.zarrurl,
        args.delete_in,
        args.sites_list,
        args.ext,
        args.chl_list,
        args.num_levels,
        args.coarsening_xy,
        args.coarsening_z,
    )
