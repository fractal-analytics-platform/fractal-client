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
    rows=None,
    cols=None,
    channels=None,
    num_levels=5,
    coarsening_xy=2,
    coarsening_z=1,
    delete_input=False,
):
    """
    Convert Yokogawa output (png, tif) to zarr file

    :param zarrurl: structure of the zarr folder
    :type zarrurl: str
    :param in_path: directory containing the input files
    :type in_path: str
    :param ext: source images extension
    :type ext: str
    :param delete_input: delete input files
    :type delete_input: bool
    :param rows: number of rows of the well
    :type rows: list
    :param cols: number of columns of the well
    :type cols: list
    :param chl_list: list of the channels   #FIXME
    :type chl_list: list                    #FIXME
    :param num_levels: number of levels in the zarr pyramid
    :type num_levels: int
    :param coarsening_xy: coarsening factor along X and Y
    :type coarsening_xy: int
    :param coarsening_z: coarsening factor along Z
    :type coarsening_z: int
    """

    if not in_path.endswith("/"):
        in_path += "/"

    # Hard-coded values (by now) of chunk sizes to be passed to rechunk,
    # both at level 0 (before coarsening) and at levels 1,2,.. (after
    # repeated coarsening).
    # Note that balance=True may override these values.
    chunk_size_x = 2560
    chunk_size_y = 2160

    # Define well
    if not zarrurl.endswith("/"):
        zarrurl += "/"
    well_row = zarrurl.split("/")[-4]
    well_column = zarrurl.split("/")[-3]
    well_ID = well_row + well_column

    lazy_imread = delayed(imread)
    fc_list = {level: [] for level in range(num_levels)}

    print(channels)

    for channel in channels:
        A, C = channel.split("_")

        l_rows = []
        all_rows = []

        print(zarrurl, well_ID)
        print(in_path + f"*_{well_ID}_*{A}_*{C}." + ext)
        filenames = sorted(
            glob(in_path + f"*_{well_ID}_*{A}*{C}." + ext), key=sort_fun
        )
        if len(filenames) == 0:
            glob_path = in_path + f"*_{well_ID}_*{A}*{C}." + ext
            raise Exception(
                "Error in yokogawa_to_zarr: len(filenames)=0.\n"
                f"  in_path: {in_path}\n"
                f"  ext: {ext}\n"
                f"  well_ID: {well_ID}\n"
                f"  channel: {channel},\n"
                f"  glob_path: {glob_path}"
            )
        max_z = max(
            [
                re.findall(r"Z(.*)C", filename.split("/")[-1])[0]
                for filename in filenames
            ]
        )

        sample = imread(filenames[0])

        start = 0
        end = int(max_z)

        for ro in range(int(rows)):
            cell = []

            for co in range(int(cols)):
                lazy_arrays = [lazy_imread(fn) for fn in filenames[start:end]]
                start += int(max_z)
                end += int(max_z)
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
        coarsening = {1: coarsening_xy, 2: coarsening_xy}
        f_matrices = {}
        for level in range(num_levels):
            if level == 0:
                f_matrices[level] = da.concatenate(all_rows, axis=1).rechunk(
                    {1: chunk_size_y, 2: chunk_size_x}, balance=True
                )
                # After concatenate, f_matrices[0] has three dimensions: z,y,x
                if coarsening_z > 1:
                    f_matrices[level] = da.coarsen(
                        np.min,
                        f_matrices[level],
                        {0: coarsening_z},
                        trim_excess=True,
                    )
            else:
                f_matrices[level] = da.coarsen(
                    np.min,
                    f_matrices[level - 1],
                    coarsening,
                    trim_excess=True,
                ).rechunk(
                    {
                        1: chunk_size_y,
                        2: chunk_size_x,
                    },
                    balance=True,
                )
            fc_list[level].append(f_matrices[level])

        if delete_input:
            for f in filenames:
                try:
                    os.remove(f)
                except OSError as e:
                    print("Error: %s : %s" % (f, e.strerror))

    level_data = [
        da.stack(fc_list[level], axis=0) for level in range(num_levels)
    ]

    shape_list = []
    for i, level in enumerate(level_data):
        level.to_zarr(zarrurl + f"{i}/", dimension_separator="/")
        print(f"Chunks at level {i}:\n", level.chunks)
        shape_list.append(level.shape)
    print()

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
        help="list of channels ",  # FIXME
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

    parser.add_argument(
        "-d",
        "--delete_input",
        action="store_true",
        help="Delete input files",
    )

    args = parser.parse_args()

    yokogawa_to_zarr(
        args.zarrurl,
        in_path=args.in_path,
        ext=args.ext,
        rows=args.rows,
        cols=args.cols,
        chl_list=args.chl_list,
        num_levels=args.num_levels,
        coarsening_xy=args.coarsening_xy,
        coarsening_z=args.coarsening_z,
        delete_input=args.delete_input,
    )
