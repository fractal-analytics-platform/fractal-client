import json

import dask.array as da

from fractal.tasks.lib_pyramid_creation import create_pyramid


def maximum_intensity_projection(
    zarrurl,
    chl_list=None,
    coarsening_xy=2,
):

    """
    Perform maximum-intensity projection along Z axis, and store the output in
    a new zarr file.

    :param zarrurl: input zarr file, at the site level (e.g. x.zarr/B/03/0/)
    :type zarrurl: str
    :param chl_list: list of channels
    :type chl_list: list
    :param coarsening_xy: coarsening factor along X and Y
    :type coarsening_z: xy

    """

    # Hard-coded values (by now) of chunk sizes to be passed to rechunk,
    # both at level 0 (before coarsening) and at levels 1,2,.. (after
    # repeated coarsening).
    # Note that balance=True may override these values.
    chunk_size_x = 1280
    chunk_size_y = 1080

    if not zarrurl.endswith("/"):
        zarrurl += "/"

    zarrurl_mip = zarrurl.replace(".zarr/", "_mip.zarr/")

    with open(zarrurl_mip + ".zattrs", "r") as inputjson:
        zattrs = json.load(inputjson)
    num_levels = len(zattrs["multiscales"][0]["datasets"])

    # Load 0-th level
    data_chl_z_y_x = da.from_zarr(zarrurl + "/0")
    # Loop over channels
    accumulate_chl = []
    for ind_chl, chl in enumerate(chl_list):
        # Perform MIP for each channel of level 0
        mip_yx = da.stack([da.max(data_chl_z_y_x[ind_chl], axis=0)], axis=0)
        accumulate_chl.append(mip_yx)
        print(chl, mip_yx.shape, mip_yx.chunks)
        print(zarrurl_mip + f"0/{ind_chl}")
    accumulate_chl = da.stack(accumulate_chl, axis=0)
    print()

    pyramid = create_pyramid(
        accumulate_chl,
        coarsening_z=1,
        coarsening_xy=coarsening_xy,
        num_levels=num_levels,
        chunk_size_x=chunk_size_x,
        chunk_size_y=chunk_size_y,
        chl_list=chl_list,
    )

    for ind_level in range(num_levels):
        pyramid[ind_level].to_zarr(
            zarrurl_mip + f"{ind_level}/", dimension_separator="/"
        )


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(prog="maximum_intensity_projection.py")
    parser.add_argument(
        "-z", "--zarrurl", help="zarr url, at the FOV level", required=True
    )

    parser.add_argument(
        "-C",
        "--chl_list",
        nargs="+",
        help="list of channels ",
    )

    parser.add_argument(
        "-cxy",
        "--coarsening_xy",
        default=2,
        type=int,
        help="coarsening factor along X and Y (optional, defaults to 2)",
    )

    args = parser.parse_args()
    maximum_intensity_projection(
        args.zarrurl, args.chl_list, args.coarsening_xy
    )
