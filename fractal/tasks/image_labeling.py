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
import time

import dask.array as da
import numpy as np
import zarr
from cellpose import core
from cellpose import models

from fractal.tasks.lib_pyramid_creation import create_pyramid_3D


def apply_label_to_single_FOV_column(
    column,
    model=None,
    do_3D=True,
    anisotropy=None,
    diameter=40.0,
    cellprob_threshold=0.0,
):

    # anisotropy = XXX  #FIXME
    mask, flows, styles, diams = model.eval(
        column,
        channels=[0, 0],
        do_3D=do_3D,
        net_avg=False,
        augment=False,
        diameter=diameter,
        anisotropy=anisotropy,
        cellprob_threshold=cellprob_threshold,
    )
    # mask = np.random.normal(10000, 200, size=column.shape).astype(np.uint16)

    return mask


def image_labeling(
    zarrurl,
    coarsening_xy=2,
    anisotropy=None,
    do_3D=True,
    diameter=None,
    cellprob_threshold=None,
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

    # Label dtype
    label_dtype = np.uint32

    # Load some level and some channel
    data_zyx = da.from_zarr(f"{zarrurl}{ref_level}")[ind_channel]

    # FIXME DO INFERENCE ABOUT 2D/3D
    do_3D = data_zyx.shape[0] > 1
    if do_3D:
        anisotropy = 5.0 / 0.325
    label_name = "label_A01_C01"
    label_axes = ["z", "y", "x"]

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

    data_zyx_rechunked = da.rechunk(
        data_zyx, chunks=(nz, img_size_y, img_size_x)
    )

    # Extract num_levels
    zattrs_file = f"{zarrurl}.zattrs"
    with open(zattrs_file, "r") as jsonfile:
        zattrs = json.load(jsonfile)
        num_levels = len(zattrs["multiscales"][0]["datasets"])
    print("num_levels", num_levels)
    print()

    # Initialize cellpose
    use_gpu = core.use_gpu()
    model = models.Cellpose(gpu=use_gpu, model_type="nuclei")

    # Initialize other things
    cumulative_shift = 0
    mask_rechunked = da.empty(
        data_zyx_rechunked.shape,
        chunks=data_zyx_rechunked.chunks,
        dtype=label_dtype,
    )
    with open("LOG_image_labeling", "w") as out:
        out.write("Total well shape/chunks:\n")
        out.write(f"{data_zyx_rechunked.shape}\n")
        out.write(f"{data_zyx_rechunked.chunks}\n\n")

    # Sequential labeling (and relabeling)
    # https://stackoverflow.com/a/72018364/19085332
    for inds in itertools.product(*map(range, mask_rechunked.blocks.shape)):

        column_data = data_zyx_rechunked.blocks[inds].compute()

        t0 = time.perf_counter()
        with open("LOG_image_labeling", "a") as out:
            out.write(f"Selected chunk: {inds}\n")
            out.write("Now running cellpose\n")
            out.write(
                f"column_data: {type(column_data)}, {column_data.shape}\n"
            )

        # This is a numpy block
        column_mask = apply_label_to_single_FOV_column(
            column_data, model=model, do_3D=do_3D, anisotropy=anisotropy
        )
        if not do_3D:
            column_mask = np.expand_dims(column_mask, axis=0)
        max_column_mask = np.max(column_mask)

        column_mask[column_mask > 0] += cumulative_shift

        cumulative_shift += max_column_mask
        if cumulative_shift > np.iinfo(label_dtype).max - 1000:
            raise Exception(
                "ERROR in re-labeling:\n"
                f"Reached more than {cumulative_shift} labels, "
                f"but dtype={label_dtype}"
            )

        t1 = time.perf_counter()
        with open("LOG_image_labeling", "a") as out:
            out.write(
                f"End, dtype={column_mask.dtype} shape={column_mask.shape}\n"
            )
            out.write(f"Elapsed: {t1-t0:.4f} seconds\n")
            out.write(
                f"column_max: {max_column_mask}, "
                f"new cumulative_shift: {cumulative_shift}\n\n"
            )

        # Put np data into a dask array.. ?? Is this out-of-memory!!  #FIXME
        start_z = inds[0] * nz
        end_z = (inds[0] + 1) * nz
        start_y = inds[1] * img_size_y
        end_y = (inds[1] + 1) * img_size_y
        start_x = inds[2] * img_size_x
        end_x = (inds[2] + 1) * img_size_x
        mask_rechunked[
            start_z:end_z, start_y:end_y, start_x:end_x
        ] = column_mask[:, :, :]

    # Rechunk to get back to the original chunking (with separate Z planes)
    mask = mask_rechunked.rechunk(data_zyx.chunks)

    # Construct resolution pyramid
    pyramid = create_pyramid_3D(
        mask,
        coarsening_z=1,
        coarsening_xy=coarsening_xy,
        num_levels=num_levels,
        chunk_size_x=img_size_x,
        chunk_size_y=img_size_y,
    )

    # Write zattrs for labels and for specific label
    # FIXME one channel? many channels? update
    labels_group = zarr.group(f"{zarrurl}labels")
    labels_group.attrs["labels"] = [label_name]
    label_group = labels_group.create_group(label_name)
    label_group.attrs["image-label"] = {"version": "0.4"}
    label_group.attrs["multiscales"] = [
        {
            "name": label_name,
            "version": "0.4",
            "axes": [
                {"name": axis_name, "type": "space"}
                for axis_name in label_axes
            ],
            "datasets": [
                {"path": f"{ind_level}"} for ind_level in range(num_levels)
            ],
        }
    ]

    # Write data into output zarr
    for ind_level in range(num_levels):
        pyramid[ind_level].astype(label_dtype).to_zarr(
            zarrurl,
            component=f"labels/{label_name}/{ind_level}",
            dimension_separator="/",
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
