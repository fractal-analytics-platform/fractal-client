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
import itertools
import json

import dask.array as da
import numpy as np
import zarr
from cellpose import core
from cellpose import models
from ome_zarr.io import parse_url
from ome_zarr.scale import Scaler
from ome_zarr.writer import write_labels


def apply_label_to_single_FOV_column(
    img,
    model=None,
):

    with open("LOG_image_labeling", "a") as out:
        out.write(f"Running cellpose on array of shape {img.shape}\n")

    # anisotropy = XXX  #FIXME

    mask, flows, styles, diams = model.eval(
        img, channels=[0, 0], do_3D=True, net_avg=False, augment=False
    )

    with open("LOG_image_labeling", "a") as out:
        out.write(f"End, dtype={mask.dtype}m shape={mask.shape}\n")

    return mask


def image_labeling(
    zarrurl,
):

    """
    FIXME
    """

    ref_level = 0  # FIXME: level for segmentation, now pinned
    ind_channel = 0  # FIXME: channel for segmentation, now pinned

    if not zarrurl.endswith("/"):
        zarrurl += "/"

    if ref_level > 0:
        raise NotImplementedError(
            "By now we can only segment the highest-resolution level"
        )

    # Load some level and some channel
    data_zyx = da.from_zarr(f"{zarrurl}{ref_level}")[ind_channel]

    # Check that input array is made of images (in terms of shape/chunks)
    img_size_y = 2160
    img_size_x = 2560
    nz, ny, nx = data_zyx.shape
    if (ny % img_size_y != 0) or (nx % img_size_x != 0):
        raise Exception(
            "Error in image_labeling, " f"data_zyx.shape: {data_zyx.shape}"
        )
    chunks_z, chunks_y, chunks_x = data_zyx.chunks
    if len(set(chunks_z)) != 1 or chunks_z[0] != 1:
        raise Exception(f"Error in image_labeling, chunks_z: {chunks_z}")
    if len(set(chunks_y)) != 1 or chunks_y[0] != img_size_y:
        raise Exception(f"Error in image_labeling, chunks_y: {chunks_y}")
    if len(set(chunks_x)) != 1 or chunks_x[0] != img_size_x:
        raise Exception(f"Error in image_labeling, chunks_x: {chunks_x}")

    use_gpu = core.use_gpu()
    model = models.Cellpose(gpu=use_gpu, model_type="nuclei")

    data_zyx_rechunked = da.rechunk(
        data_zyx, chunks=(nz, img_size_y, img_size_x)
    )
    mask_rechunked = data_zyx_rechunked.map_blocks(
        apply_label_to_single_FOV_column,
        chunks=(nz, img_size_y, img_size_x),
        meta=np.array((), dtype=np.uint16),
        model=model,
    )

    # Relabeling via explicit loop
    # (https://stackoverflow.com/a/72018364/19085332)
    cumulative_shift = 0
    for inds in itertools.product(*map(range, mask_rechunked.blocks.shape)):
        max_label = da.max(mask_rechunked.blocks[inds]).compute()
        print("max_label, cumulative_shift", max_label, cumulative_shift)
        start_z = inds[0] * nz
        end_z = (inds[0] + 1) * nz
        start_y = inds[1] * img_size_y
        end_y = (inds[1] + 1) * img_size_y
        start_x = inds[2] * img_size_x
        end_x = (inds[2] + 1) * img_size_x
        mask_rechunked[
            start_z:end_z, start_y:end_y, start_x:end_x
        ] += cumulative_shift
        cumulative_shift += max_label

    # Rechunk
    mask = da.rechunk(mask_rechunked, chunks=(1, img_size_y, img_size_x))

    # From here on, it's all about storing the mask

    # Extract num_levels
    zattrs_file = f"{zarrurl}.zattrs"
    with open(zattrs_file, "r") as jsonfile:
        zattrs = json.load(jsonfile)
        num_levels = len(zattrs["multiscales"][0]["datasets"])
    print("num_levels", num_levels)
    print()

    cxy = 2  # Still needed, for the scaler
    # FIXME: check that shapes change by x2 each time

    # Define custom ome_zarr scaler
    our_scaler = Scaler()
    our_scaler.downscale = cxy
    our_scaler.max_layer = num_levels - 1

    # FIXME: relabeling

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
        scaler=our_scaler,
    )


def old_image_labeling(
    zarrurl,
):

    """
    FIXME
    """

    raise NotImplementedError(
        "This is an old function, which we only keep here because it includes"
        " the logic for rescaling when segmenting a level>0"
    )

    ref_level = 3  # FIXME: level for segmentation, now pinned
    ind_channel = 0  # FIXME: channel for segmentation, now pinned

    # if ref_level > 0:

    # Load some level and some channel
    data_zyx = da.from_zarr(f"{zarrurl}/{ref_level}")[ind_channel]
    shape_zyx_segmented_level = data_zyx.shape

    with open("LOG_image_labeling", "w") as out:
        out.write(f"shape of data to segment: {shape_zyx_segmented_level}\n")
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
