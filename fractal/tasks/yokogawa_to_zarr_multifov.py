import os
import re
from glob import glob

import dask.array as da
import numpy as np
from dask import delayed
from skimage.io import imread

from fractal.tasks.lib_parse_filename_metadata import parse_metadata


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
    num_levels=5,
    coarsening_xy=2,
    coarsening_z=1,
    delete_input=False,
):

    """
    Convert Yokogawa output (png, tif) to zarr file, keeping the multi-site
    (AKA multi-field-of-view) structure in the zarr file.

    """

    if not in_path.endswith("/"):
        in_path += "/"
    if not zarrurl.endswith("/"):
        zarrurl += "/"

    # WARNING: the zarr file should point to a well, not to a plate
    if zarrurl.endswith(".zarr/"):
        raise Exception(
            "Error in yokogawa_to_zarr_multifov, "
            f"zarrurl={zarrurl} should be at the site level."
        )

    # Define well and site
    if not zarrurl.endswith("/"):
        zarrurl += "/"
    well_row = zarrurl.split("/")[-4]
    well_column = zarrurl.split("/")[-3]
    well_ID = well_row + well_column
    site = zarrurl.split("/")[-2]
    print(f"zarrurl={zarrurl}")
    print(f"site={site}")
    if not site.isdigit():
        raise Exception(
            "Error in yokogawa_to_zarr_multifov\n"
            f"zarrurl={zarrurl}\n"
            f"site={site}"
        )
    site = int(site)
    site += 1  # FIXME: is this robust?

    lazy_imread = delayed(imread)

    print(f"Channels: {chl_list}")

    matrix_level_channel = {}

    for ind_chl, chl in enumerate(chl_list):
        A, C = chl.split("_")

        print(zarrurl, well_ID)
        # Find all filenames for a given well and cannel, then filter by site
        # [note: this is to avoid specifying the format of site (e.g. "F001")
        glob_path = f"{in_path}*_{well_ID}_*{A}*{C}.{ext}"
        filenames = sorted(
            [
                fn
                for fn in glob(glob_path)
                if int(parse_metadata(fn.split("/")[-1])["F"]) == site
            ]
        )
        # FIXME: should we also specify plate?

        if len(filenames) == 0:
            raise Exception(
                "Error in yokogawa_to_zarr: len(filenames)=0.\n"
                f"  in_path: {in_path}\n"
                f"  ext: {ext}\n"
                f"  well_ID: {well_ID}\n"
                f"  site: {site}\n"
                f"  channel: {chl},\n"
                f"  glob_path: {glob_path}"
            )
        print(f"[SITE {site}, CHANNEL {chl}] {len(filenames)} matching images")

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
                        np.mean,
                        matrices_site_channel_levels[level],
                        {0: coarsening_z},
                        trim_excess=True,
                    )
            else:
                matrices_site_channel_levels.append(
                    da.coarsen(
                        np.mean,
                        matrices_site_channel_levels[level - 1],
                        coarsening,
                        trim_excess=True,
                    )
                )

            matrix_level_channel[
                (level, ind_chl)
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
        level_stack.to_zarr(zarrurl + f"{level}/", dimension_separator="/")

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
        num_levels=args.num_levels,
        coarsening_xy=args.coarsening_xy,
        coarsening_z=args.coarsening_z,
        delete_input=args.delete_input,
    )
