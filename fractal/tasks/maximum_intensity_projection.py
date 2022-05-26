import json

import dask.array as da
import numpy as np


def maximum_intensity_projection(
    zarrurl,
    chl_list=None,
    coarsening_factor_xy=2,
):

    """
    #FIXME docstring
    # zarrul = xxx.zarr/B/03/0/
    """

    # Hard-coded values (by now) of chunk sizes to be passed to rechunk,
    # both at level 0 (before coarsening) and at levels 1,2,.. (after
    # repeated coarsening).
    # Note that balance=True may override these values.
    chunk_size_x = 256 * 5
    chunk_size_y = 216 * 5

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

    # Create coarser levels
    coarsening = {1: coarsening_factor_xy, 2: coarsening_factor_xy}
    all_channels_per_level = {level: [] for level in range(num_levels)}
    for ind_chl, chl in enumerate(chl_list):
        f_matrices = {}
        for level in range(num_levels):
            if level == 0:
                # f_matrices[level] shape: zyx (with dummy z)
                f_matrices[level] = accumulate_chl[ind_chl].rechunk(
                    {1: chunk_size_x, 2: chunk_size_y}, balance=True
                )
            else:
                f_matrices[level] = da.coarsen(
                    np.min, f_matrices[level - 1], coarsening, trim_excess=True
                ).rechunk(
                    {1: max(1, chunk_size_x), 2: max(1, chunk_size_y)},
                    balance=True,
                )

            all_channels_per_level[level].append(f_matrices[level])

    # all_channels_per_level[level] shape: czyx
    for ind_level in range(num_levels):
        level = da.stack(all_channels_per_level[ind_level], axis=0)
        level.to_zarr(zarrurl_mip + f"{ind_level}/", dimension_separator="/")
    print()


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
