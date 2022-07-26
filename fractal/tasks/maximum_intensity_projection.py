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
import json

import dask.array as da

from fractal.tasks.lib_pyramid_creation import write_pyramid


def maximum_intensity_projection(
    zarrurl,
    coarsening_xy=2,
):

    """
    Perform maximum-intensity projection along Z axis, and store the output in
    a new zarr file.

    :param zarrurl: input zarr file, at the site level (e.g. x.zarr/B/03/0/)
    :type zarrurl: str
    :param coarsening_xy: coarsening factor along X and Y
    :type coarsening_xy: xy

    """

    # Hard-coded values (by now) of chunk sizes to be passed to rechunk,
    # both at level 0 (before coarsening) and at levels 1,2,.. (after
    # repeated coarsening).
    # Note that balance=True may override these values.
    img_size_x = 2560
    img_size_y = 2160

    if not zarrurl.endswith("/"):
        zarrurl += "/"

    zarrurl_mip = zarrurl.replace(".zarr/", "_mip.zarr/")

    with open(zarrurl_mip + ".zattrs", "r") as inputjson:
        zattrs = json.load(inputjson)
    num_levels = len(zattrs["multiscales"][0]["datasets"])

    # Load 0-th level
    data_chl_z_y_x = da.from_zarr(zarrurl + "/0")
    num_channels = data_chl_z_y_x.shape[0]
    # Loop over channels
    accumulate_chl = []
    for ind_ch in range(num_channels):

        # Perform MIP for each channel of level 0
        mip_yx = da.stack([da.max(data_chl_z_y_x[ind_ch], axis=0)], axis=0)

        accumulate_chl.append(mip_yx)
        print(ind_ch, mip_yx.shape, mip_yx.chunks)
        print(zarrurl_mip + f"0/{ind_ch}")
    accumulate_chl = da.stack(accumulate_chl, axis=0)
    print()

    # Construct resolution pyramid
    write_pyramid(
        accumulate_chl,
        newzarrurl=zarrurl_mip,
        overwrite=False,
        coarsening_xy=coarsening_xy,
        num_levels=num_levels,
        chunk_size_x=img_size_x,
        chunk_size_y=img_size_y,
    )


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(prog="maximum_intensity_projection.py")
    parser.add_argument(
        "-z", "--zarrurl", help="zarr url, at the FOV level", required=True
    )

    parser.add_argument(
        "-cxy",
        "--coarsening_xy",
        default=2,
        type=int,
        help="coarsening factor along X and Y (optional, defaults to 2)",
    )

    args = parser.parse_args()
    maximum_intensity_projection(args.zarrurl, args.coarsening_xy)
