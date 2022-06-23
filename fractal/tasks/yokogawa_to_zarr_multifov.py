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


def yokogawa_to_zarr_multifov(
    zarrurl,
    in_path=None,
    ext=None,
    chl_list=None,
    sites_dict={},
    num_levels=5,
    coarsening_xy=2,
    coarsening_z=1,
    delete_input=False,
):

    """
    Convert Yokogawa output (png, tif) to zarr file, keeping the multi-site
    (AKA multi-field-of-view) structure in the zarr file.

    """

    sites_list = sites_dict[zarrurl]

    if not in_path.endswith("/"):
        in_path += "/"
    if not zarrurl.endswith("/"):
        zarrurl += "/"

    # WARNING: the zarr file should point to a well, not to a plate
    if zarrurl.endswith(".zarr/"):
        raise Exception(
            "Error in replicate_zarr_structure, "
            f"zarrurl={zarrurl} does not end with .zarr/"
        )

    r = zarrurl.split("/")[-3]
    c = zarrurl.split("/")[-2]

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
                1: coarsening_xy,  # Y dimension
                2: coarsening_xy,  # X dimension
            }
            for level in range(num_levels):
                if level == 0:
                    matrices_site_channel_levels = [matrix_site_channel]
                    if coarsening_z > 1:
                        matrices_site_channel_levels[level] = da.coarsen(
                            np.min,
                            matrices_site_channel_levels[level],
                            {0: coarsening_z},
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

        if delete_input:
            for f in filenames:
                try:
                    os.remove(f)
                except OSError as e:
                    print("Error: %s : %s" % (f, e.strerror))

        for level in range(num_levels):
            level_list = [
                matrix_level_channel[key]
                for key in matrix_level_channel.keys()
                if key[0] == level
            ]
            level_stack = da.stack(level_list, axis=0)
            level_stack.to_zarr(
                zarrurl_site + f"{level}/", dimension_separator="/"
            )

    # FIXME: return something useful
    shape_list = []

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
        help="structure of the zarr folder (at the well level)",
    )

    parser.add_argument(
        "-d",
        "--delete_input",
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

    yokogawa_to_zarr_multifov(
        args.zarrurl,
        in_path=args.in_path,
        ext=args.ext,
        chl_list=args.chl_list,
        sites_list=args.sites_list,
        num_levels=args.num_levels,
        coarsening_xy=args.coarsening_xy,
        coarsening_z=args.coarsening_z,
        delete_input=args.delete_input,
    )
