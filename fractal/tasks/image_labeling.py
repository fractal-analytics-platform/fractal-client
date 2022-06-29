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
import dask.array as da
import zarr
from cellpose import core
from cellpose import models
from ome_zarr.io import parse_url
from ome_zarr.writer import write_labels


def image_labeling(
    zarrurl,
):

    """
    FIXME
    """

    # Load some level and some channel
    # FIXME: this is now pinned to level 0 and channel 0
    dapi_dset = da.from_zarr(zarrurl + "/0")[0]

    # Some debugging info
    print("LOADED ZARR")
    print(dapi_dset.shape)
    print()

    # Perform segmentation
    use_gpu = core.use_gpu()
    model = models.Cellpose(gpu=use_gpu, model_type="nuclei")
    mask, flows, styles, diams = model.eval(
        dapi_dset, channels=[0, 0], do_3D=True, net_avg=False, augment=False
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
    label_axes = [
        "z",
        "y",
        "x",
    ]  # could change if e.g. a 2D image was processed => can we infer this?
    write_labels(mask, group=root, name=label_name, axes=label_axes)


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
