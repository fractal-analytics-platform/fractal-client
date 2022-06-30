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
import zarr
from cellpose import core
from cellpose import models
from ome_zarr.io import parse_url
from ome_zarr.scale import Scaler
from ome_zarr.writer import write_labels


def image_labeling(
    zarrurl,
):

    """
    FIXME
    """

    ref_level = 2  # FIXME: level for segmentation, now pinned
    ind_channel = 0  # FIXME: channel for segmentation, now pinned

    # Load some level and some channel
    data_zyx = da.from_zarr(f"{zarrurl}/{ref_level}")[ind_channel]
    shape_zyx_segmented_level = data_zyx.shape

    # Some debugging info
    print("LOADED ZARR", shape_zyx_segmented_level)
    print()

    # Extract num_levels
    zattrs_file = f"{zarrurl}.zattrs"
    with open(zattrs_file, "r") as jsonfile:
        zattrs = json.load(jsonfile)
        num_levels = len(zattrs["multiscales"][0]["datasets"])
    print("num_levels", num_levels)
    print()

    cxy = 2  # Still needed, for the scaler
    # FIXME: check that shapes change by x2 each time

    # Define custom ome_zarr coordinate transformations
    coordinatetransformations = []
    for level in range(num_levels):

        zarray_file = f"{zarrurl}{level}/.zarray"
        with open(zarray_file, "r") as jsonfile:
            shape_zxy = json.load(jsonfile)["shape"][1:]
        if level == 0:
            shape_zxy_highres = shape_zxy[:]
            global_prefactor_y = (
                shape_zxy_highres[1] / shape_zyx_segmented_level[1]
            )
            global_prefactor_x = (
                shape_zxy_highres[2] / shape_zyx_segmented_level[2]
            )

        factor_y = global_prefactor_y * shape_zxy_highres[1] / shape_zxy[1]
        factor_x = global_prefactor_x * shape_zxy_highres[2] / shape_zxy[2]
        factors = [1.0, factor_y, factor_x]

        # Alternative way
        # factor = cxy ** (level + ref_level) * 1.0
        # factors = [1.0, factor, factor]

        print(f"level {level}, factors {factors}")
        print()

        coordinatetransformations.append(
            [
                {
                    "type": "scale",
                    "scale": factors,
                }
            ]
        )

    # Define custom ome_zarr scaler
    our_scaler = Scaler()
    our_scaler.downscale = cxy
    our_scaler.max_layer = num_levels - 1

    # Perform segmentation
    use_gpu = core.use_gpu()
    model = models.Cellpose(gpu=use_gpu, model_type="nuclei")
    mask, flows, styles, diams = model.eval(
        data_zyx, channels=[0, 0], do_3D=True, net_avg=False, augment=False
    )

    # Some debugging info
    print("CELLPOSE OUTPUT")
    print(type(mask), mask.shape)
    print(type(flows), len(flows), len(flows[0]))
    print(type(styles), styles.shape)
    print(type(diams), diams)
    print()

    # Store labels
    store = parse_url(zarrurl, mode="w").store
    root = zarr.group(store=store)
    label_name = "label_image"
    label_axes = ["z", "y", "x"]  # FIXME: should be inferred (2D vs 3D)

    write_labels(
        mask,
        group=root,
        name=label_name,
        axes=label_axes,
        coordinate_transformations=coordinatetransformations,
        scaler=our_scaler,
    )


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(prog="image_labeling.py")
    parser.add_argument(
        "-z", "--zarrurl", help="zarr url, at the FOV level", required=True
    )

    args = parser.parse_args()
    image_labeling(
        args.zarrurl,
    )
