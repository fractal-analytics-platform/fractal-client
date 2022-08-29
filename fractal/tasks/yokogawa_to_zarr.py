"""
Copyright 2022 (C) Friedrich Miescher Institute for Biomedical Research and
University of Zurich

Original authors:
Tommaso Comparin <tommaso.comparin@exact-lab.it>
Marco Franzon <marco.franzon@exact-lab.it>

This file is part of Fractal and was originally developed by eXact lab S.r.l.
<exact-lab.it> under contract with Liberali Lab from the Friedrich Miescher
Institute for Biomedical Research and Pelkmans Lab from the University of
Zurich.
"""
import os
import re
from glob import glob
from pathlib import Path
from typing import Any
from typing import Dict
from typing import Iterable
from typing import Optional

import dask.array as da
from dask import delayed
from skimage.io import imread

from fractal.tasks.lib_pyramid_creation import write_pyramid


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
    *,
    input_paths: Iterable[Path],
    output_path: Path,
    rows: int,
    cols: int,
    delete_input=False,
    metadata: Optional[Dict[str, Any]] = None,
    component: str = None,
):
    """
    Convert Yokogawa output (png, tif) to zarr file

    Example arguments:
      input_paths[0] = /tmp/input/*png  (Path)
      output_path = /tmp/output/*zarr   (Path)
      metadata = {"channel_list": [...], "num_levels": ..., }
      component = plate.zarr/B/03/0/
    """

    from devtools import debug

    debug(output_path)
    debug(component)

    chl_list = metadata["channel_list"]
    original_path_list = metadata["original_paths"]
    in_path = Path(original_path_list[0]).parent
    ext = Path(original_path_list[0]).name
    num_levels = metadata["num_levels"]
    coarsening_xy = metadata["coarsening_xy"]

    # Hard-coded values (by now) of chunk sizes to be passed to rechunk,
    # both at level 0 (before coarsening) and at levels 1,2,.. (after
    # repeated coarsening).
    # Note that balance=True may override these values.
    chunk_size_x = 2560
    chunk_size_y = 2160

    # Define well
    component_split = component.split("/")
    well_row = component_split[1]
    well_column = component_split[2]

    well_ID = well_row + well_column

    delayed_imread = delayed(imread)

    print(f"Channels: {chl_list}")

    list_channels = []
    for chl in chl_list:
        A, C = chl.split("_")

        glob_path = f"{in_path}/*_{well_ID}_*{A}*{C}{ext}"
        print(f"glob path: {glob_path}")
        filenames = sorted(glob(glob_path), key=sort_fun)
        if len(filenames) == 0:
            raise Exception(
                "Error in yokogawa_to_zarr: len(filenames)=0.\n"
                f"  in_path: {in_path}\n"
                f"  ext: {ext}\n"
                f"  well_ID: {well_ID}\n"
                f"  channel: {chl},\n"
                f"  glob_path: {glob_path}"
            )
        max_z = max(
            [
                re.findall(r"Z(.*)C", filename.split("/")[-1])[0]
                for filename in filenames
            ]
        )

        sample = imread(filenames[0])

        # Loop over FOVs, corresponding to filenames[start:end]
        start = 0
        end = int(max_z)
        data_zfyx = []
        for row in range(int(rows)):
            FOV_rows = []
            for col in range(int(cols)):
                images = [delayed_imread(fn) for fn in filenames[start:end]]
                lazy_images = [
                    da.from_delayed(
                        image, shape=sample.shape, dtype=sample.dtype
                    )
                    for image in images
                ]
                z_stack = da.stack(lazy_images, axis=0)
                start += int(max_z)
                end += int(max_z)
                FOV_rows.append(z_stack)
            data_zfyx.append(da.block(FOV_rows))
        # Remove FOV index
        data_zyx = da.concatenate(data_zfyx, axis=1)
        list_channels.append(data_zyx)
    data_czyx = da.stack(list_channels, axis=0)

    if delete_input:
        for f in filenames:
            try:
                os.remove(f)
            except OSError as e:
                print("Error: %s : %s" % (f, e.strerror))

    # Construct resolution pyramid
    write_pyramid(
        data_czyx,
        newzarrurl=output_path.parent.as_posix() + f"/{component}",
        overwrite=False,
        coarsening_xy=coarsening_xy,
        num_levels=num_levels,
        chunk_size_x=chunk_size_x,
        chunk_size_y=chunk_size_y,
    )


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
        help="list of channel names (e.g. A01_C01)",
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
        delete_input=args.delete_input,
    )
