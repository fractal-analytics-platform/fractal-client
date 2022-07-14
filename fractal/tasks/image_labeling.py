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
import time
from concurrent.futures import ThreadPoolExecutor

import dask
import dask.array as da
import numpy as np
from cellpose import core
from cellpose import models

from fractal.tasks.lib_pyramid_creation import write_pyramid


def segment_FOV(
    column,
    block_info=None,
    model=None,
    do_3D=True,
    anisotropy=None,
    diameter=40.0,
    cellprob_threshold=0.0,
    label_dtype=None,
):

    chunk_location = block_info[None]["chunk-location"]

    # Write some debugging info
    with open("LOG_image_labeling", "a") as out:
        out.write(
            f"[{chunk_location}] START Cellpose |"
            f" column: {type(column)}, {column.shape} |"
            f" do_3D: {do_3D}\n"
        )

    # Actual labeling
    t0 = time.perf_counter()
    mask, flows, styles, diams = model.eval(
        column,
        channels=[0, 0],
        do_3D=do_3D,
        net_avg=False,
        augment=False,
        diameter=diameter,
        anisotropy=anisotropy,
        cellprob_threshold=cellprob_threshold,
        # batch_size=16,
    )
    if not do_3D:
        mask = np.expand_dims(mask, axis=0)
    t1 = time.perf_counter()

    # Write some debugging info
    with open("LOG_image_labeling", "a") as out:
        out.write(
            f"[{chunk_location}] END   Cellpose |"
            f" Elapsed: {t1-t0:.4f} seconds |"
            f" mask shape: {mask.shape},"
            f" mask dtype: {mask.dtype} (before recast to {label_dtype}),"
            f" max(mask): {np.max(mask)}\n"
        )
        # out.write(
        #    f"[{chunk_location}] "
        #    f" flows shape: {np.array(flows[1]).shape},"
        #    f" flows dtype: {np.array(flows[1]).dtype}\n")

    return mask.astype(label_dtype)


def image_labeling(
    zarrurl,
    coarsening_xy=2,
    labeling_level=0,
    labeling_channel=None,
    chl_list=None,
    num_threads=1,
    relabeling=False,
    # More parameters
    anisotropy=None,
    diameter=None,
    cellprob_threshold=None,
):

    """
    FIXME
    """

    # Sanitize zarr path
    if not zarrurl.endswith("/"):
        zarrurl += "/"

    # Find channel index
    if labeling_channel not in chl_list:
        raise Exception(f"ERROR: {labeling_channel} not in {chl_list}")
    ind_channel = chl_list.index(labeling_channel)

    # Check that level=0
    if labeling_level > 0:
        raise NotImplementedError(
            "By now we can only segment the highest-resolution level"
        )

    # Set labels dtype
    label_dtype = np.uint32

    # Load ZYX data
    data_zyx = da.from_zarr(f"{zarrurl}{labeling_level}")[ind_channel]

    # Select 2D/3D behavior and set some parameters
    do_3D = data_zyx.shape[0] > 1
    if do_3D:
        if anisotropy is None:
            # Reasonable value for level 0 (for some of our UZH datasets)
            anisotropy = 1.0 / 0.1625

    # Load .zattrs file
    zattrs_file = f"{zarrurl}.zattrs"
    with open(zattrs_file, "r") as jsonfile:
        zattrs = json.load(jsonfile)

    # Extract num_levels
    num_levels = len(zattrs["multiscales"][0]["datasets"])
    print("num_levels", num_levels)
    print()

    # Try to read channel label from OMERO metadata
    try:
        omero_label = zattrs["omero"]["channels"][ind_channel]["label"]
        label_name = f"label_{omero_label}"
    except (KeyError, IndexError):
        label_name = f"label_{ind_channel}"

    # Check that input array is made of images (in terms of shape/chunks)
    img_size_y = 2160
    img_size_x = 2560
    nz, ny, nx = data_zyx.shape
    if (ny % img_size_y != 0) or (nx % img_size_x != 0):
        raise Exception(
            "Error in image_labeling, data_zyx.shape: {data_zyx.shape}"
        )
    chunks_z, chunks_y, chunks_x = data_zyx.chunks
    if len(set(chunks_z)) != 1 or chunks_z[0] != 1:
        raise Exception(f"Error in image_labeling, chunks_z: {chunks_z}")
    if len(set(chunks_y)) != 1 or chunks_y[0] != img_size_y:
        raise Exception(f"Error in image_labeling, chunks_y: {chunks_y}")
    if len(set(chunks_x)) != 1 or chunks_x[0] != img_size_x:
        raise Exception(f"Error in image_labeling, chunks_x: {chunks_x}")

    # Initialize cellpose
    use_gpu = core.use_gpu()
    model = models.Cellpose(gpu=use_gpu, model_type="nuclei")

    # Initialize other things
    with open("LOG_image_labeling", "w") as out:
        out.write(f"Start image_labeling task for {zarrurl}\n")
        out.write(f"relabeling: {relabeling}\n")
        out.write("Total well shape/chunks:\n")
        out.write(f"{data_zyx.shape}\n")
        out.write(f"{data_zyx.chunks}\n\n")

    # Map labeling function onto all chunks (i.e., FOV colums)
    mask = (
        data_zyx.rechunk((nz, img_size_y, img_size_x))
        .map_blocks(
            segment_FOV,
            meta=np.array((), dtype=label_dtype),
            model=model,
            do_3D=do_3D,
            anisotropy=anisotropy,
            label_dtype=label_dtype,
        )
        .rechunk((1, img_size_y, img_size_x))
    )

    with open("LOG_image_labeling", "a") as out:
        out.write(
            f"After map_block, mask will have shape {mask.shape} "
            f"and chunks {mask.chunks}\n\n"
        )

    # Construct resolution pyramid
    with dask.config.set(pool=ThreadPoolExecutor(num_threads)):
        write_pyramid(
            mask,
            newzarrurl=f"{zarrurl}labels/{label_name}/",
            overwrite=False,
            coarsening_xy=coarsening_xy,
            num_levels=num_levels,
            chunk_size_x=2560,
            chunk_size_y=2160,
        )


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(prog="image_labeling.py")
    parser.add_argument(
        "-z", "--zarrurl", help="zarr url, at the FOV level", required=True
    )
    parser.add_argument(
        "-C",
        "--chl_list",
        nargs="+",
        help="list of channel names (e.g. A01_C01)",
    )
    parser.add_argument(
        "-cxy",
        "--coarsening_xy",
        default=2,
        type=int,
        help="coarsening factor along X and Y (optional, defaults to 2)",
    )
    parser.add_argument(
        "-lc",
        "--labeling_channel",
        help="name of channel for labeling (e.g. A01_C01)",
    )
    parser.add_argument(
        "--num_threads",
        default=1,
        type=int,
        help="TBD",
    )

    args = parser.parse_args()
    image_labeling(
        args.zarrurl,
        coarsening_xy=args.coarsening_xy,
        chl_list=args.chl_list,
        labeling_channel=args.labeling_channel,
        num_threads=args.num_threads,
        # FIXME: more arguments
    )
