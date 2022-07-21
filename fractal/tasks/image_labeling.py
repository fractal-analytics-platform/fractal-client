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
import shutil
import time
from concurrent.futures import ThreadPoolExecutor

import anndata as ad
import dask
import dask.array as da
import numpy as np
import zarr
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
    logfile="LOG_image_labeling",
):

    chunk_location = block_info[None]["chunk-location"]

    # Write some debugging info
    with open(logfile, "a") as out:
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
    )
    if not do_3D:
        mask = np.expand_dims(mask, axis=0)
    t1 = time.perf_counter()

    # Write some debugging info
    with open(logfile, "a") as out:
        out.write(
            f"[{chunk_location}] END   Cellpose |"
            f" Elapsed: {t1-t0:.4f} seconds |"
            f" mask shape: {mask.shape},"
            f" mask dtype: {mask.dtype} (before recast to {label_dtype}),"
            f" max(mask): {np.max(mask)}\n"
        )

    return mask.astype(label_dtype)


def image_labeling(
    zarrurl,
    coarsening_xy=2,
    labeling_level=0,
    labeling_channel=None,
    chl_list=None,
    num_threads=1,
    relabeling=True,
    anisotropy=None,
    diameter=None,
    cellprob_threshold=None,
    model_type="nuclei",
):

    """
    FIXME
    """

    # Sanitize zarr path
    if not zarrurl.endswith("/"):
        zarrurl += "/"

    # Find well ID
    # FIXME: only useful for our temporary log files
    well_ID = "_".join(zarrurl.split("/")[-4:-2])
    logfile = f"LOG_image_labeling_{well_ID}"

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

            adata = ad.read_zarr(f"{zarrurl}tables/FOV_ROI_table")
            pixel_size_x = adata["pixel_size_x"][0]
            pixel_size_y = adata["pixel_size_y"][0]
            pixel_size_z = adata["pixel_size_z"][0]
            if not np.allclose(pixel_size_x, pixel_size_y):
                raise Exception(
                    "ERROR: XY anisotropy detected\n"
                    f"pixel_size_x={pixel_size_x}\n"
                    f"pixel_size_y={pixel_size_y}"
                )
            anisotropy = pixel_size_z / pixel_size_x
            if labeling_level > 0:
                raise NotImplementedError(
                    "TODO: fix automatic anisotropy "
                    "detection for higher levels"
                )

    # Check model_type
    if model_type not in ["nuclei", "cyto2", "cyto"]:
        raise Exception(f"ERROR model_type={model_type} is not allowed.")

    # Load zattrs file
    zattrs_file = f"{zarrurl}.zattrs"
    with open(zattrs_file, "r") as jsonfile:
        zattrs = json.load(jsonfile)

    # Preliminary checks on multiscales
    multiscales = zattrs["multiscales"]
    if len(multiscales) > 1:
        raise Exception(f"ERROR: There are {len(multiscales)} multiscales")
    if "coordinateTransformations" in multiscales[0].keys():
        raise Exception(
            "ERROR: coordinateTransformations at the multiscales "
            "level are not currently supported"
        )

    # Extract num_levels
    num_levels = len(multiscales[0]["datasets"])
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
    model = models.Cellpose(gpu=use_gpu, model_type=model_type)

    # Initialize other things
    with open(logfile, "w") as out:
        out.write(f"Start image_labeling task for {zarrurl}\n")
        out.write(f"relabeling: {relabeling}\n")
        out.write(f"model_type: {model_type}\n")
        out.write(f"num_threads: {num_threads}\n")
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
            diameter=diameter,
            cellprob_threshold=cellprob_threshold,
            logfile=logfile,
        )
        .rechunk((1, img_size_y, img_size_x))
    )

    with open(logfile, "a") as out:
        out.write(
            f"After map_block, mask will have shape {mask.shape} "
            f"and chunks {mask.chunks}\n\n"
        )

    # Write zattrs for labels and for specific label
    # FIXME deal with: (1) many channels, (2) overwriting
    labels_group = zarr.group(f"{zarrurl}labels")
    labels_group.attrs["labels"] = [label_name]
    label_group = labels_group.create_group(label_name)
    label_group.attrs["image-label"] = {"version": "0.4"}
    label_group.attrs["multiscales"] = [
        {
            "name": label_name,
            "version": "0.4",
            "axes": multiscales[0]["axes"],
            "datasets": multiscales[0]["datasets"],
        }
    ]

    if relabeling:

        # Execute all, and write level-0 mask to disk
        with dask.config.set(pool=ThreadPoolExecutor(num_threads)):
            write_pyramid(
                mask,
                newzarrurl=f"{zarrurl}labels/{label_name}/",
                overwrite=False,
                coarsening_xy=coarsening_xy,
                num_levels=1,
                chunk_size_x=2560,
                chunk_size_y=2160,
            )

        with open(logfile, "a") as out:
            out.write("\nStart relabeling\n")
        t0 = time.perf_counter()

        # Load non-relabeled mask from disk
        mask = da.from_zarr(zarrurl, component=f"labels/{label_name}/{0}")
        mask_rechunked = mask.rechunk((nz, img_size_y, img_size_x))
        newmask_rechunked = da.empty(
            shape=mask_rechunked.shape,
            chunks=mask_rechunked.chunks,
            dtype=label_dtype,
        )

        # Sequential relabeling
        # https://stackoverflow.com/a/72018364/19085332
        num_labels_tot = 0
        num_labels_column = 0
        for inds in itertools.product(
            *map(range, mask_rechunked.blocks.shape)
        ):

            # Select a specific chunk (=column in 3D, =image in 2D)
            column_mask = mask_rechunked.blocks[inds].compute()
            num_labels_column = np.max(column_mask)

            # Apply re-labeling and update total number of labels
            shift = np.zeros_like(column_mask, dtype=label_dtype)
            shift[column_mask > 0] = num_labels_tot
            num_labels_tot += num_labels_column

            with open(logfile, "a") as out:
                out.write(
                    f"Chunk {inds}, "
                    f"num_labels_column={num_labels_column}, "
                    f"num_labels_tot={num_labels_tot}\n"
                )

            # Check that total number of labels is under control
            if num_labels_tot > np.iinfo(label_dtype).max:
                raise Exception(
                    "ERROR in re-labeling:\n"
                    f"Reached {num_labels_tot} labels, "
                    f"but dtype={label_dtype}"
                )
            # Re-assign the chunk to the new array
            start_z = inds[0] * nz
            end_z = (inds[0] + 1) * nz
            start_y = inds[1] * img_size_y
            end_y = (inds[1] + 1) * img_size_y
            start_x = inds[2] * img_size_x
            end_x = (inds[2] + 1) * img_size_x
            newmask_rechunked[start_z:end_z, start_y:end_y, start_x:end_x] = (
                column_mask + shift
            )

        newmask = newmask_rechunked.rechunk(data_zyx.chunks)

        # FIXME: this is ugly
        shutil.rmtree(zarrurl + f"labels/{label_name}/{0}")

        write_pyramid(
            newmask,
            newzarrurl=f"{zarrurl}labels/{label_name}/",
            overwrite=False,
            coarsening_xy=coarsening_xy,
            num_levels=num_levels,
            chunk_size_x=2560,
            chunk_size_y=2160,
        )

        t1 = time.perf_counter()
        with open(logfile, "a") as out:
            out.write(f"End relabeling, elapsed: {t1-t0} s\n")

    else:

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
                aggregation_function=np.max,
            )

        with open(logfile, "a") as out:
            out.write("\nSkip relabeling\n")


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
