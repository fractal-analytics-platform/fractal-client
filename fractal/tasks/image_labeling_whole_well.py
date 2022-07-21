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

import dask.array as da
import numpy as np
import zarr
from cellpose import core
from cellpose import models

from fractal.tasks.lib_pyramid_creation import write_pyramid


def image_labeling_whole_well(
    zarrurl,
    coarsening_xy=2,
    labeling_level=2,
    labeling_channel=None,
    chl_list=None,
    relabeling=True,
    anisotropy=None,
    diameter_level0=35.0,
    cellprob_threshold=0.0,
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
    logfile = f"LOG_image_labeling_{well_ID}_whole_well"

    # Check model_type
    if model_type not in ["nuclei", "cyto2", "cyto"]:
        raise Exception(f"ERROR model_type={model_type} is not allowed.")

    # Work on MIP zarr file
    # FIXME: this is a temporary hack
    zarrurl = zarrurl.replace(".zarr/", "_mip.zarr/")

    # Find channel index
    if labeling_channel not in chl_list:
        raise Exception(f"ERROR: {labeling_channel} not in {chl_list}")
    ind_channel = chl_list.index(labeling_channel)

    # Check that level>=1
    if labeling_level < 1:
        raise NotImplementedError("By now we can only segment levels >= 1")

    # Load ZYX data
    data_zyx = da.from_zarr(f"{zarrurl}{labeling_level}")[ind_channel]

    # Select 2D/3D behavior and set some parameters
    if data_zyx.shape[0] > 1:
        raise Exception(
            f"ERROR shape = {data_zyx.shape}"
            " but there can be only one Z plane."
        )

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

    # Extract axes, and remove channel
    new_axes = [ax for ax in multiscales[0]["axes"] if ax["type"] != "channel"]

    # Try to read channel label from OMERO metadata
    try:
        omero_label = zattrs["omero"]["channels"][ind_channel]["label"]
        label_name = f"label_{omero_label}"
    except (KeyError, IndexError):
        label_name = f"label_{ind_channel}"

    # Initialize cellpose
    use_gpu = core.use_gpu()
    model = models.Cellpose(gpu=use_gpu, model_type=model_type)

    # Initialize other things
    with open(logfile, "w") as out:
        out.write(f"Start image_labeling_whole_well task for {zarrurl}\n")
        out.write(f"model_type: {model_type}\n")
        out.write("Total well shape/chunks:\n")
        out.write(f"{data_zyx.shape}\n")
        out.write(f"{data_zyx.chunks}\n\n")

    # Write some debugging info
    with open(logfile, "a") as out:
        out.write(
            f"START Cellpose |" f" well: {type(data_zyx)}, {data_zyx.shape} \n"
        )

    # Actual labeling
    t0 = time.perf_counter()
    mask, flows, styles, diams = model.eval(
        data_zyx,
        channels=[0, 0],
        do_3D=False,
        net_avg=False,
        augment=False,
        diameter=(diameter_level0 / coarsening_xy**labeling_level),
        cellprob_threshold=cellprob_threshold,
    )
    mask = np.expand_dims(mask, axis=0)
    t1 = time.perf_counter()

    # Write some debugging info
    with open(logfile, "a") as out:
        out.write(
            f"END   Cellpose |"
            f" Elapsed: {t1-t0:.4f} seconds |"
            f" mask shape: {mask.shape},"
            f" mask dtype: {mask.dtype},"
            f" max(mask): {np.max(mask)}\n\n"
        )

    # Convert mask to dask
    mask_da = da.from_array(mask).rechunk((1, 2160, 2560))
    with open(logfile, "a") as out:
        out.write(
            f"da.from_array(upscaled_mask) [with rechunking]: {mask_da}\n\n"
        )

    # Construct rescaled datasets
    datasets = multiscales[0]["datasets"]
    new_datasets = []
    for ds in datasets:
        new_ds = {}

        # Copy all keys that are not coordinateTransformations (e.g. path)
        for key in ds.keys():
            if key != "coordinateTransformations":
                new_ds[key] = ds[key]

        # Update coordinateTransformations
        old_transformations = ds["coordinateTransformations"]
        new_transformations = []
        for t in old_transformations:
            if t["type"] == "scale":
                new_t = {"type": "scale"}
                new_t["scale"] = [
                    t["scale"][0],
                    t["scale"][1] * coarsening_xy**labeling_level,
                    t["scale"][2] * coarsening_xy**labeling_level,
                ]
                new_transformations.append(new_t)
            else:
                new_transformations.append(t)
        new_ds["coordinateTransformations"] = new_transformations

        new_datasets.append(new_ds)

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
            "axes": new_axes,
            "datasets": new_datasets,
        }
    ]

    # Construct resolution pyramid
    write_pyramid(
        mask_da,
        newzarrurl=f"{zarrurl}labels/{label_name}/",
        overwrite=False,
        coarsening_xy=coarsening_xy,
        num_levels=num_levels,
        chunk_size_x=2560,
        chunk_size_y=2160,
        aggregation_function=np.max,
    )


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(prog="image_labeling_whole_well.py")
    parser.add_argument(
        "-z",
        "--zarrurl",
        help="zarr url, at the merged-FOV level",
        required=True,
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

    args = parser.parse_args()
    image_labeling_whole_well(
        args.zarrurl,
        coarsening_xy=args.coarsening_xy,
        chl_list=args.chl_list,
        labeling_channel=args.labeling_channel,
    )
