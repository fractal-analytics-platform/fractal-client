import dask.array as da


def maximum_intensity_projection(
    zarrurl,
    chl_list=None,
):

    """
    #FIXME docstring
    # zarrul = xxx.zarr/B/03/0/
    """

    if not zarrurl.endswith("/"):
        zarrurl += "/"

    zarrurl_mip = zarrurl.replace(".zarr/", "_mip.zarr/")

    # Load 0-th level
    data_chl_z_y_x = da.from_zarr(zarrurl + "/0")

    # Loop over channels
    accumulate_chl = []
    for ind_chl, chl in enumerate(chl_list):

        mip_yx = da.stack([da.max(data_chl_z_y_x[ind_chl], axis=0)], axis=0)
        accumulate_chl.append(mip_yx)
        print(chl, mip_yx.shape, mip_yx.chunks)
        print(zarrurl_mip + f"0/{ind_chl}")
    accumulate_chl = da.stack(accumulate_chl, axis=0)
    accumulate_chl.to_zarr(zarrurl_mip + "0/", dimension_separator="/")


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

    args = parser.parse_args()
    maximum_intensity_projection(args.zarrurl, args.chl_list)


"""
    z_list = [rc.split("/")[-1] for rc in glob(zarrurl + "*/*/*")]
    y_list = [rc.split("/")[-1] for rc in glob(zarrurl + "*/*/*/*")]
    x_list = [rc.split("/")[-1] for rc in glob(zarrurl + "*/*/*/*/*")]
    z_unique = sorted(list(set(z_list)))
    y_unique = sorted(list(set(y_list)))
    x_unique = sorted(list(set(x_list)))
data_zyx = []
for z in z_unique:
    data_yx = []
    for y in y_unique:
        data_x = []
        for x in x_unique:
            tmp_data_x = da.from_zarr(...)
            data_x.append(tmp_data_x)
            #lazy_data_x = lazy_load(zarrurl + f"0/{int(chl) - 1}/{z}/{y}/{x}")
            #data_x.append(da.from_delayed(lazy_data_x, dtype="<f8"))   #
        data_yx.append(da.concatenate(data_x))
    data_zyx.append(da.stack(data_yx, axis=0))
data_zyx = da.stack(data_zyx, axis=0)

#FIXME: perform maximization at each z level, rather than at the end
"""
